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

from optimade.filterparser import LarkParser
from optimade.filtertransformers.elasticsearch import Transformer, Quantity


class FilterException(Exception):
    ''' Raised on parsing a filter expression with syntactic of semantic errors. '''
    pass


_quantities: Dict[str, Quantity] = None
_parser = LarkParser(version=(0, 10, 1))
_transformer = None


def parse_filter(filter_str: str) -> Q:
    ''' Parses the given optimade filter str and returns a suitable elastic search query.

    Arguments:
        filter_str: Can be direct user input with no prior processing.

    Raises:
        FilterException: If the given str cannot be parsed, or if there are any semantic
            errors in the given expression.
    '''
    global _quantities
    global _transformer
    if _quantities is None:
        from nomad.datamodel import OptimadeEntry
        _quantities = {
            q.name: Quantity(
                q.name, es_field='dft.optimade.%s' % q.name,
                elastic_mapping_type=q.a_search.mapping.__class__)

            for q in OptimadeEntry.m_def.all_quantities.values()
            if 'search' in q.m_annotations}

        _quantities['elements'].length_quantity = _quantities['nelements']
        _quantities['dimension_types'].length_quantity = _quantities['dimension_types']
        _quantities['elements'].has_only_quantity = Quantity(name='only_atoms')
        _quantities['elements'].nested_quantity = _quantities['elements_ratios']
        _quantities['elements_ratios'].nested_quantity = _quantities['elements_ratios']

        _transformer = Transformer(quantities=_quantities.values())

    try:
        parse_tree = _parser.parse(filter_str)
    except Exception as e:
        raise FilterException('Syntax error: %s' % str(e))

    try:
        query = _transformer.transform(parse_tree)
    except Exception as e:
        raise FilterException('Semantic error: %s' % str(e))

    return query
