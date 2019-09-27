# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from optimade.filterparser import LarkParser
import lark
from elasticsearch_dsl import Q, Text, Keyword, Integer
import ase.data

from nomad.metainfo.optimade import OptimadeStructureEntry
from nomad.metainfo import Quantity


class FilterException(Exception):
    """ Raised on parsing a filter expression with syntactic of semantic errors. """
    pass


_cmp_operators = {'>': 'gt', '>=': 'gte', '<': 'lt', '<=': 'lte'}
_rev_cmp_operators = {'>': '<', '>=': '<=', '<': '>', '<=': '=>'}
_has_operators = {'ALL': 'must', 'ANY': 'should'}
_length_quantities = {'elements': 'nelements', 'elements_rations': 'nelements', 'dimension_types': 'dimension_types'}


class Transformer(lark.Transformer):
    """ Transformer for the Lark parser generator used for the filterparser.

    It translates the parse tree into an elastic search query.
    """

    def _field(self, quantity, nested=None):
        optimade_field_name = quantity.name
        if nested is not None:
            optimade_field_name = '%s.%s' % (nested, optimade_field_name)
        return 'optimade.%s' % optimade_field_name

    def _order_terms(self, l, o, r):
        if isinstance(l, Quantity):
            if isinstance(r, Quantity):
                raise Exception('Cannot compare two quantities: %s, %s' % (l.name, r.name))

            return l, o, r
        else:
            if isinstance(r, Quantity):
                o = _rev_cmp_operators.get(o, o)
                return r, o, l

            raise Exception('Cannot compare two values: %s, %s' % (str(l), str(l)))

    def _query(self, quantity, o, value, nested=None):
        field = self._field(quantity, nested=nested)
        if o in _cmp_operators:
            return Q('range', **{field: {_cmp_operators[o]: value}})

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

        return Q('nested', path='optimade.elements_ratios', query=dict(bool=dict(must=queries)))

    def _wildcard_query(self, quantity, wildcard):
        return Q('wildcard', **{self._field(quantity): wildcard})

    def __default__(self, tree, children, *args, **kwargs):
        """ Default behavior for rules that only replace one symbol with another """
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
        o, = args
        return ~o

    def cmp_op(self, args):
        l, o, r = args
        field, o, value = self._order_terms(l, o, r)
        return self._query(field, o, value)

    def has_op(self, args):
        quantities, predicates = args
        return self._has_query(quantities, predicates)

    def has_list_op(self, args):
        quantities, o, predicates_list = args
        queries = [
            self._has_query(quantities, predicates)
            for predicates in predicates_list]

        if o in _has_operators:
            return Q('bool', **{_has_operators[o]: queries})

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
        if quantity.name not in _length_quantities:
            raise Exception('LENGTH is not supported for %s' % quantity.name)

        return OptimadeStructureEntry.m_section.quantities[_length_quantities[quantity.name]]

    def known_op(self, args):
        quantity, qualifier = args
        query = Q('exists', field=self._field(quantity))
        if qualifier == 'KNOWN':
            return query
        elif qualifier == 'UNKNOWN':
            return ~query  # pylint: disable=invalid-unary-operand-type

        raise NotImplementedError

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


_parser = LarkParser(version=(0, 10, 0))
_transformer = Transformer()


def parse_filter(filter_str: str) -> Q:
    """ Parses the given optimade filter str and returns a suitable elastic search query.

    Arguments:
        filter_str: Can be direct user input with no prior processing.

    Raises:
        FilterException: If the given str cannot be parsed, or if there are any semantic
            errors in the given expression.
    """

    try:
        parse_tree = _parser.parse(filter_str)
    except Exception as e:
        raise FilterException('Syntax error: %s' % str(e))

    try:
        query = _transformer.transform(parse_tree)
    except Exception as e:
        raise FilterException('Semantic error: %s' % str(e))

    return query
