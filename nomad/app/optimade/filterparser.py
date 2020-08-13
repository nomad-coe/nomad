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

from typing import Dict
from elasticsearch_dsl import Q
from cachetools import cached

from optimade.filterparser import LarkParser
from optimade.filtertransformers.elasticsearch import ElasticTransformer, Quantity

from nomad.search import search_quantities


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

    quantities['elements'].length_quantity = quantities['nelements']
    quantities['dimension_types'].length_quantity = quantities['dimension_types']
    quantities['elements'].has_only_quantity = Quantity(name='only_atoms')
    quantities['elements'].nested_quantity = quantities['elements_ratios']
    quantities['elements_ratios'].nested_quantity = quantities['elements_ratios']

    if nomad_properties is not None:
        for search_quantity in search_quantities.values():
            name = search_quantity.name
            if '.' in name:
                if name.startswith(nomad_properties):
                    name = name[len(nomad_properties) + 1:]
                else:
                    continue

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
