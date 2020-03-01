from typing import Callable, Any

from nomad import metainfo


# TODO multi, split are more flask related
class SearchQuantity:
    '''
    A metainfo quantity annotation class that defines additional properties that determine
    how to search for the respective quantity. Only quantities that have this will
    be mapped to elastic search.

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
        default_statistic: Indicates this quantity to be part of the default statistics.
        statistics_size:
            The maximum number of values in a statistic. Default is 10.
        group: Indicates that his quantity can be used to group results. The value will
            be the name of the group.
        es_quantity: The quantity in the elastic mapping that is used to search. This is
            especially useful if the quantity represents a inner document and only one
            quantity of this inner object is used. Default is the name of the quantity.
        es_mapping: A valid elasticsearch_dsl mapping. Default is ``Keyword()``.
        es_value: A callable that is applied to section to get a value for this quantity in the elastic index.
        derived: A callable that is applied to search parameter values before search.
    '''

    def __init__(
            self,
            name: str = None, description: str = None,
            many_and: str = None, many_or: str = None,
            order_default: bool = False,
            group: str = None, metric: str = None, metric_name: str = None,
            default_statistic: bool = False,
            statistic_size: int = 10,
            es_quantity: str = None,
            es_mapping: Any = None,
            es_value: Callable[[Any], Any] = None,
            derived: Callable[[Any], Any] = None):

        self.name = name
        self.description = description
        self.many_and = many_and
        self.many_or = many_or
        self.order_default = order_default
        self.group = group
        self.default_statistic = default_statistic
        self.metric = metric
        self.metric_name = metric_name
        self.statistic_size = statistic_size
        self.es_quantity = es_quantity
        self.es_mapping = es_mapping
        self.es_value = es_value
        self.derived = derived

        self.prefix: str = None
        self.qualified_name: str = None

        assert many_and is None or many_or is None, 'A search quantity can only be used for multi or many search'
        assert many_and in [None, 'split', 'append'], 'Only split and append are valid values'
        assert many_or in [None, 'split', 'append'], 'Only split and append are valid values'

    def configure(self, quantity: metainfo.Quantity, prefix: str = None):
        if self.name is None:
            self.name = quantity.name

        if self.description is None:
            self.description = quantity.description

        if prefix is not None:
            self.qualified_name = '%s.%s' % (prefix, self.name)
            if self.es_quantity is not None:
                self.es_quantity = '%s.%s' % (prefix, self.es_quantity)
            if self.metric_name is not None:
                self.metric_name = '%s.%s' % (prefix, self.metric_name)
            if self.group is not None:
                self.group = '%s.%s' % (prefix, self.group)
        else:
            self.qualified_name = self.name

        if self.es_quantity is None:
            self.es_quantity = self.qualified_name
        if self.metric_name is None and self.metric is not None:
            self.metric_name = self.qualified_name

    @property
    def argparse_action(self):
        if self.many_or is not None:
            return self.many_or

        if self.many_and is not None:
            return self.many_and

        return None

    @property
    def many(self):
        return self.many_and is not None or self.many_or is not None


def init(section: metainfo.MSection):
    pass
