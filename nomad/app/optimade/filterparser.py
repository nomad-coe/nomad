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
from optimade.filterparser import LarkParser
from optimade.filtertransformers.elasticsearch import Transformer, Quantity
from elasticsearch_dsl import Q
from nomad.metainfo.optimade import OptimadeEntry


class FilterException(Exception):
    """ Raised on parsing a filter expression with syntactic of semantic errors. """
    pass


quantities: Dict[str, Quantity] = {
    q.name: Quantity(
        q.name, es_field='optimade.%s' % q.name,
        elastic_mapping_type=q.m_annotations['elastic']['type'])

    for q in OptimadeEntry.m_def.all_quantities.values()
    if 'elastic' in q.m_annotations}

quantities['elements'].length_quantity = quantities['nelements']
quantities['dimension_types'].length_quantity = quantities['dimension_types']
quantities['elements'].has_only_quantity = Quantity(name='only_atoms')
quantities['elements'].nested_quantity = quantities['elements_ratios']
quantities['elements_ratios'].nested_quantity = quantities['elements_ratios']


_parser = LarkParser(version=(0, 10, 0))
_transformer = Transformer(quantities=quantities.values())


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
