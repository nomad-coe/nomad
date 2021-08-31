#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from typing import Dict
from elasticsearch_dsl import Q, Date
from cachetools import cached

from optimade.filterparser import LarkParser
from optimade.filtertransformers.elasticsearch import (
    Quantity, ElasticTransformer as OPTElasticTransformer)
from optimade.models import CHEMICAL_SYMBOLS, ATOMIC_NUMBERS

from .common import provider_specific_fields


_parser = LarkParser(version=(0, 10, 1))


class FilterException(Exception):
    ''' Raised on parsing a filter expression with syntactic of semantic errors. '''
    pass


@cached(cache={})
def _get_transformer(nomad_properties, without_prefix):
    from nomad.datamodel import OptimadeEntry
    quantities: Dict[str, Quantity] = {
        q.name: Quantity(
            q.name, es_field='dft.optimade.%s' % q.name,
            elastic_mapping_type=q.a_search.mapping.__class__)

        for q in OptimadeEntry.m_def.all_quantities.values()
        if 'search' in q.m_annotations}

    quantities['id'] = Quantity('id', es_field='calc_id')
    quantities['immutable_id'] = Quantity('immutable_id', es_field='calc_id')
    quantities['last_modified'] = Quantity(
        'last_modified', es_field='upload_time', elastic_mapping_type=Date)

    quantities['elements'].length_quantity = quantities['nelements']
    quantities['elements'].has_only_quantity = Quantity(name='only_atoms')
    quantities['elements'].nested_quantity = quantities['elements_ratios']
    quantities['elements_ratios'].nested_quantity = quantities['elements_ratios']

    if nomad_properties is not None:
        for name, search_quantity in provider_specific_fields():
            names = ['_nmd_' + name]
            if without_prefix:
                names.append(name)

            for name in names:
                if name not in quantities:
                    quantities[name] = Quantity(
                        name,
                        es_field=search_quantity.search_field,
                        elastic_mapping_type=search_quantity.mapping.__class__)

    return ElasticTransformer(quantities=quantities.values())


def parse_filter(filter_str: str, nomad_properties='dft', without_prefix=False) -> Q:
    ''' Parses the given optimade filter str and returns a suitable elastic search query.

    Arguments:
        filter_str: Can be direct user input with no prior processing.
        nomad_properties: Also include the nomad proprietary properties of the given domain.
        without_prefix: Do not prefix the nomad proprietary properties with _nmd_.

    Raises:
        FilterException: If the given str cannot be parsed, or if there are any semantic
            errors in the given expression.
    '''
    transformer = _get_transformer(nomad_properties, without_prefix)

    try:
        parse_tree = _parser.parse(filter_str)
    except Exception as e:
        raise FilterException('Syntax error: %s' % str(e))

    try:
        query = transformer.transform(parse_tree)
    except Exception as e:
        raise FilterException('Semantic error: %s' % str(e))

    return query


class ElasticTransformer(OPTElasticTransformer):
    def _has_query_op(self, quantities, op, predicate_zip_list):
        # We override this to add 'HAS ONLY' support.
        if op == 'HAS ONLY':
            # HAS ONLY comes with heavy limitations, because there is no such thing
            # in elastic search. Only supported for elements, where we can construct
            # an anonymous 'formula' based on elements sorted by order number and
            # can do a = comparision to check if all elements are contained
            if len(quantities) > 1:
                raise Exception('HAS ONLY is not supported with zip')
            quantity = quantities[0]

            if quantity.has_only_quantity is None:
                raise Exception('HAS ONLY is not supported by %s' % quantity.name)

            def values():
                for predicates in predicate_zip_list:
                    if len(predicates) != 1:
                        raise Exception('Tuples not supported in HAS ONLY')
                    op, value = predicates[0]
                    if op != '=':
                        raise Exception('Predicated not supported in HAS ONLY')
                    if not isinstance(value, str):
                        raise Exception('Only strings supported in HAS ONLY')
                    yield value

            try:
                order_numbers = list([ATOMIC_NUMBERS[element] for element in values()])
                order_numbers.sort()
                value = ''.join(
                    [CHEMICAL_SYMBOLS[number - 1] for number in order_numbers]
                )
            except KeyError:
                raise NotImplementedError('HAS ONLY is only supported for chemical symbols')

            return Q('term', **{quantity.has_only_quantity.name: value})

        else:
            return super()._has_query_op(quantities, op, predicate_zip_list)

    def property_zip_addon(self, args):
        return args

    def value_zip(self, args):
        return self.value_list(args)

    def value_zip_list(self, args):
        return args
