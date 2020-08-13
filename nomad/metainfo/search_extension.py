# Copyright 2020 Markus Scheidgen
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

from typing import Callable, Any, Dict, List

from nomad.metainfo.elastic_extension import Elastic


search_quantities: Dict[str, 'Search'] = {}
''' All available search quantities by their full qualified name. '''

metrics: Dict[str, 'Search'] = {}
'''
The available search metrics. Metrics are integer values given for each entry that can
be used in statistics (aggregations), e.g. the sum of all total energy calculations or
cardinality of all unique geometries. The key is the metric name.
'''

groups: Dict[str, 'Search'] = {}
''' The available groupable quantities. The key is the group name. '''

order_default_quantities: Dict[str, 'Search'] = {}
''' The quantity for each domain (key) that is the default quantity to order search results by. '''


# TODO multi, split are more flask related
class Search(Elastic):
    '''
    A metainfo quantity annotation class that defines additional properties that determine
    how to search for the respective quantity. Only quantities that have this will
    be mapped to elastic search. The annotation is an extension of :class:`Elastic` and
    add nomad API specific search features like grouping, statistics, metrics, domains, etc.

    Attributes:
        name: The name of this search quantity. Will be the name in the elastic index and
            the name for the search parameter. Default is the metainfo quantity name.
        many_or: Indicates that an 'or' (es terms) search is performed if many values are given.
            Otherwise an 'and' (es bool->should->match) is performed.  Values are 'split' and
            'append' to indicate how URL search parameters should be treated.
        many_and: Indicates that many values can be supplied for search. Values are 'split' and
            'append' to indicate how URL search parameters should be treated.
        order_default: Indicates that this quantity is used to order search results
            if no other ordering was specificed.
        metric: Quantity can be used to build statistics. Statistics provide a metric
            value for each value of the quantity. E.g. number of datasets with a given atom label.
            This defines a metric based on this quantity. Values need to be a valid
            elastic search aggregation (e.g. sum, cardinality, etc.).
        metric_name: If this quantity is indicated to function as a metric, the metric
            needs a name. By default the quantities name is used.
        statistics_size:
            The maximum number of values in a statistic. Default is 10.
        statistics_order:
            The order key that is passed to elastic search to determine the order of
            the statistic values.
        statistics_values:
            If the statistics has a fixed set of values, use this parameter to set it
            as a list of strings. Will fill statistics_size with the len of this list.
            The information can be used (e.g. by the GUI) to fill in empty values.
        group: Indicates that his quantity can be used to group results. The value will
            be the name of the group.
        derived: A callable that is applied to search parameter values before search.
        search_field: The qualified field in the elastic mapping that is used to search.
            This might be different from the field that is used to store the value in
            elastic search. This is especially useful if the field represents a inner
            document and a subfield of this inner object should be used for search.
    '''

    def __init__(
            self,
            name: str = None, description: str = None,
            many_and: str = None, many_or: str = None,
            order_default: bool = False,
            group: str = None, metric: str = None, metric_name: str = None,
            statistic_size: int = None,
            statistic_order: str = '_key',
            statistic_values: List[str] = None,
            derived: Callable[[Any], Any] = None,
            search_field: str = None,
            **kwargs):

        super().__init__(field=None, **kwargs)

        self.name = name
        self.description = description
        self.many_and = many_and
        self.many_or = many_or
        self.order_default = order_default
        self.group = group
        self.metric = metric
        self.metric_name = metric_name

        self.statistic_fixed_size = statistic_size
        self.statistic_size = statistic_size if statistic_size is not None else 20
        self.statistic_order = statistic_order
        self.statistic_values = statistic_values
        self.search_field = search_field

        self.derived = derived

        self.qualified_name: str = None

        assert many_and is None or many_or is None, 'A search quantity can only be used for multi or many search'
        assert many_and in [None, 'split', 'append'], 'Only split and append are valid values'
        assert many_or in [None, 'split', 'append'], 'Only split and append are valid values'

    def init_annotation(self, definition):
        if self.name is None:
            self.name = definition.name
        assert self.name is not None

        if self.description is None:
            self.description = definition.description

        super().init_annotation(definition)

    def register(self, prefix, field):
        domain_or_all = self.definition.m_parent.m_get_annotations('domain', '__all__')

        prefix_and_dot = prefix + '.' if prefix is not None else ''

        self.qualified_name = prefix_and_dot + self.name
        if self.search_field is not None:
            self.search_field = prefix_and_dot + self.search_field
        else:
            self.search_field = self.qualified_name

        assert self.qualified_name not in search_quantities, 'Search quantities must have a unique name: %s' % self.name
        search_quantities[self.qualified_name] = self

        if self.metric is not None:
            if self.metric_name is None:
                self.metric_name = self.qualified_name
            else:
                self.metric_name = prefix_and_dot + self.metric_name

            assert self.metric_name not in metrics, 'Metric names must be unique: %s' % self.metric_name
            metrics[self.metric_name] = self

        if self.group is not None:
            self.group = prefix_and_dot + self.group
            assert self.group not in groups, 'Groups must be unique'
            groups[self.group] = self

        if self.order_default:
            assert order_default_quantities.get(domain_or_all) is None, 'Only one quantity can be the order default'
            order_default_quantities[domain_or_all] = self

    @property
    def argparse_action(self):
        if self.many_or is not None:
            return self.many_or

        if self.many_and is not None:
            return self.many_and

        return None

    @property
    def many(self) -> bool:
        return self.many_and is not None or self.many_or is not None

    @property
    def flask_field(self):
        from flask_restplus import fields
        value_field = fields.String
        if self.definition.type == int:
            value_field = fields.Integer

        if self.many:
            return fields.List(value_field(), description=self.description)
        else:
            return value_field(description=self.description)

    @property
    def statistic_values(self):
        return self._statistic_values

    @statistic_values.setter
    def statistic_values(self, values):
        self._statistic_values = values
        if self._statistic_values is not None:
            self.statistic_size = len(self._statistic_values)
            self.statistic_fixed_size = len(self._statistic_values)
