import json
from optimade.filterparser import LarkParser
from lark import Transformer
from elasticsearch_dsl import Q, Text, Keyword, Integer
import ase.data

from nomad.processing import Upload
from nomad.search import SearchRequest
from nomad.metainfo.optimade import OptimadeStructureEntry
from nomad.metainfo import Quantity


def test_get_entry(published: Upload):
    calc_id = list(published.calcs)[0].calc_id

    with published.upload_files.archive_file(calc_id) as f:
        data = json.load(f)

    assert 'OptimadeStructureEntry' in data
    search_result = SearchRequest().search_parameter('calc_id', calc_id).execute_paginated()['results'][0]
    assert 'optimade' in search_result


class ESTransformer(Transformer):

    cmp_operators = {'>': 'gt', '>=': 'gte', '<': 'lt', '<=': 'lte'}
    has_operators = {'ALL': 'must', 'ANY': 'should'}
    length_quantities = {'elements': 'nelements', 'elements_rations': 'nelements', 'dimension_types': 'dimension_types'}

    def _field(self, quantity, nested=None):
        optimade_field_name = quantity
        if nested is not None:
            optimade_field_name = '%s.%s' % (nested, optimade_field_name)
        return 'optimade.%s' % optimade_field_name

    def _order_terms(self, l, r):
        if isinstance(l, Quantity):
            if isinstance(r, Quantity):
                raise Exception('Cannot compare two quantities: %s, %s' % (l.name, r.name))
            else:
                return l, r
        else:
            if isinstance(r, Quantity):
                return r, l
            else:
                raise Exception('Cannot compare two values: %s, %s' % (str(l), str(l)))

    def __default__(self, tree, children, *args, **kwargs):
        return children[0]

    def and_expr(self, args):
        if len(args) == 1:
            return args[0]
        l, r = args
        return l & r

    def or_expr(self, args):
        if len(args) == 1:
            return args[0]
        l, r = args
        return l | r

    def not_expr(self, args):
        if len(args) == 1:
            return args[0]
        o, = args
        return ~o

    def _query(self, quantity, o, value, nested=None):
        field = self._field(quantity, nested=nested)
        if o in ESTransformer.cmp_operators:
            return Q('range', **{field: {ESTransformer.cmp_operators[o]: value}})

        elastic_annotation = quantity.m_annotations.get('elastic', None)
        if elastic_annotation['type'] == Text:
            query_type = 'match'
        elif elastic_annotation['type'] in [Keyword, Integer]:
            query_type = 'term'
        else:
            raise NotImplementedError('Quantity has unsupported ES field type')

        if o in ['=', '']:
            return Q(query_type, **{field: value})

        if o == '!=':
            return ~Q(query_type, **{field: value})  # pylint: disable=invalid-unary-operand-type

        raise Exception('Unknown operator %s' % o)

    def cmp_op(self, args):
        l, o, r = args
        field, value = self._order_terms(l, r)
        return self._query(field, o, value)

    def has_op(self, args):
        quantities, predicates = args
        return self._has_query(quantities, predicates)

    def _has_query(self, quantities, predicates):
        if len(quantities) != len(predicates):
            raise Exception(
                'Tuple length does not match: %s <o> %s ' %
                (':'.join(quantities), ':'.join(predicates)))

        if len(quantities) == 1:
            o, value = predicates[0]
            return self._query(quantities[0], o, value)

        if any(quantity.name not in ['elements', 'elements_ratios'] for quantity in quantities):
            raise Exception('Expression with tuples are only supported for elements and elements_positions')

        queries = [
            self._query(field, o, value, nested='elements_ratios')
            for field, (o, value) in zip(quantities, predicates)]

        return Q('nested', path='elements_ratios', query=dict(bool=dict(must=queries)))

    def has_list_op(self, args):
        quantities, o, predicates_list = args
        queries = [
            self._has_query(quantities, predicates)
            for predicates in predicates_list]

        if o in ESTransformer.has_operators:
            return Q('bool', **{ESTransformer.has_operators[o]: queries})

        raise Exception('Unknown operator %s' % o)

    def has_only_op(self, args):
        quantity, lst = args

        if quantity.name != 'elements':
            raise Exception('HAS ONLY is only supported for elements')

        def values():
            for predicates in lst:
                if len(predicates) != 1:
                    raise Exception('Tuples not supported in HAS ONLY')
                op, value = predicates[0]
                if op != '':
                    raise Exception('Predicated not supported in HAS ONLY')
                if not isinstance(value, str):
                    raise Exception('Only strings supported in HAS ONLY')
                yield value

        try:
            order_numbers = list([ase.data.atomic_numbers[element] for element in values()])
            order_numbers.sort()
            value = ''.join([ase.data.chemical_symbols[number] for number in order_numbers])
        except KeyError as e:
            raise Exception('Not a chemical symbol: %s' % str(e))

        return Q('term', only_atoms=value)

    def length(self, args):
        quantity, = args
        if quantity.name not in ESTransformer.length_quantities:
            raise Exception('LENGTH is not supported for %s' % quantity.name)

        return OptimadeStructureEntry.m_section.quantities[ESTransformer.length_quantities[quantity.name]]

    def known_op(self, args):
        quantity, qualifier = args
        query = Q('exists', field=self._field(quantity))
        if qualifier == 'KNOWN':
            return query
        elif qualifier == 'UNKNOWN':
            return ~query  # pylint: disable=invalid-unary-operand-type

        raise NotImplementedError

    def _wildcard_query(self, quantity, wildcard):
        return Q('wildcard', **{self._field(quantity): dict(value=wildcard)})

    def contains_op(self, args):
        quantity, value = args
        return self._wildcard_query(quantity, '*%s*' % value)

    def starts_op(self, args):
        quantity, value = args
        return self._wildcard_query(quantity, '%s*' % value)

    def ends_op(self, args):
        quantity, value = args
        return self._wildcard_query(quantity, '*%s' % value)

    def list(self, args):
        return list(args)

    def quantity_tuple(self, args):
        return list(args)

    def predicate_tuple(self, args):
        return list(args)

    def predicate(self, args):
        if len(args) == 1:
            return '', args[0]
        else:
            return args[0], args[1]

    def quantity(self, args):
        quantity_name = args[0]
        quantity_def = OptimadeStructureEntry.m_section.quantities.get(quantity_name, None)

        if quantity_def is None:
            raise Exception('%s is not a known quantity' % quantity_name)

        elastic_annotation = quantity_def.m_annotations.get('elastic', None)
        if elastic_annotation is None:
            raise Exception('%s is not supported in queries' % quantity_name)

        return quantity_def

    def int_literal(self, args):
        return int(args[0])

    def float_literal(self, args):
        return float(args[0])

    def string_literal(self, args):
        return args[0].strip('"')


def test_optimade_parser(published: Upload):
    p = LarkParser(version=(0, 10, 0))
    tree = p.parse('''
        LENGTH elements > 2 AND
        elements:elements_ratios HAS ALL "H":>0.66,"H":<0.67 AND
        elements:elements_ratios:elements_ratios HAS ALL "O":>0.33:<0.34 AND
        (chemical_formula_reduced IS UNKNOWN OR chemical_formula_reduced CONTAINS "H2") AND
        elements HAS ONLY "O", "H" AND
        LENGTH dimension_types = 0
    ''')
    transformer = ESTransformer()
    query = transformer.transform(tree)

    result = SearchRequest(query=query).execute_paginated()
    print(result)
