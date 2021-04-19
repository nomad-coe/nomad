from typing import Callable
from lark import v_args
from elasticsearch_dsl import Q, Field, Text, Keyword, Integer, Boolean
from optimade.filtertransformers.elasticsearch import Quantity
from nomad.atomutils import get_hill_decomposition

from nomad.app.optimade.filterparser import ElasticTransformer


_cmp_operators = {">": "gt", ">=": "gte", "<": "lt", "<=": "lte"}
_rev_cmp_operators = {">": "<", ">=": "<=", "<": ">", "<=": "=>"}
_has_operators = {"ALL": "must", "ANY": "should"}
_length_quantities = {
    "elements": "nelements",
    "elements_rations": "nelements",
    "dimension_types": "dimension_types",
}


class MQuantity(Quantity):
    """A specialized version of
    optimade.filtertransformers.elasticsearch.Quantity for handling material
    queries. Added support for nested properties and converter functions.

    Attributes:
        name: The name of the quantity as used in the filter expressions.
        es_field: The name of the field for this quanity in elastic search, will be
            ``name`` by default.
        elastic_mapping_type: A decendent of an elasticsearch_dsl Field that denotes which
            mapping was used in the elastic search index.
        length_quantity: Elasticsearch does not support length of arrays, but we can
            map fields with array to other fields with ints about the array length. The
            LENGTH operator will only be supported for quantities with this attribute.
        has_only_quantity: Elasticsearch does not support exclusive search on arrays, like
            a list of chemical elements. But, we can order all elements by atomic number
            and use a keyword field with all elements to perform this search. This only
            works for elements (i.e. labels in ``CHEMICAL_SYMBOLS``) and quantities
            with this attribute.
        nested_quantity: To support optimade's 'zipped tuple' feature (e.g.
            'elements:elements_ratios HAS "H":>0.33), we use elasticsearch nested objects
            and nested queries. This quantity will provide the field for the nested
            object that contains the quantity (and others). The zipped tuples will only
            work for quantities that share the same nested object quantity.
        converter: Optional function that is used for manipulating searched
            values. Takes as input the token and returns the value that will be fed
            into the actual search.
        nested_path: The path used for nested queries. Provide only for nested
            fields.
    """

    def __init__(
        self,
        name,
        es_field: str = None,
        elastic_mapping_type: Field = None,
        length_quantity: "MQuantity" = None,
        has_only_quantity: "MQuantity" = None,
        nested_quantity: "MQuantity" = None,
        converter: Callable = None,
        nested_path: str = None,
    ):
        super().__init__(name, es_field, elastic_mapping_type, length_quantity, has_only_quantity, nested_quantity)
        self.converter = converter
        self.nested_path = nested_path


class MElasticTransformer(ElasticTransformer):
    '''
    A specialized Optimade/Lark transformer for handling material queries.
    Provides mostly the same functionality as
    optimade.filtertransformers.elasticsearch.ElasticTransformer, but has
    additions that make nested queries, parameter conversions and the use of
    boolean values possible.

    Uses elasticsearch_dsl and will produce a :class:`Q` instance.

    Arguments:
        quantities: A list of :class:`MQuantity`s that describe how optimade (and other)
            quantities are mapped to the elasticsearch index.
    '''
    def _query_op(self, quantity, op, value, nested=None):
        """
        Return a range, match, or term query for the given quantity, comparison
        operator, and value
        """
        field = self._field(quantity, nested=nested)
        if op in _cmp_operators:
            return Q("range", **{field: {_cmp_operators[op]: value}})

        if quantity.elastic_mapping_type == Text:
            query_type = "match"
        elif quantity.elastic_mapping_type in [Keyword, Integer, Boolean]:
            query_type = "term"
        else:
            raise NotImplementedError("Quantity has unsupported ES field type")

        if op in ["=", ""]:
            return Q(query_type, **{field: value})

        if op == "!=":
            return ~Q(  # pylint: disable=invalid-unary-operand-type
                query_type, **{field: value}
            )

    def _has_query_op(self, quantities, op, predicate_zip_list):
        """
        Returns a bool query that combines the operator queries ():func:`_query_op`)
        for each predicate and zipped quantity pericates combinations.
        """
        if op == "HAS":
            kind = "must"  # in case of HAS we do a must over the "list" of the one given element
        elif op == "HAS ALL":
            kind = "must"
        elif op == "HAS ANY":
            kind = "should"
        elif op == "HAS ONLY":
            # HAS ONLY comes with heavy limitations, because there is no such thing
            # in elastic search. Only supported for elements, where we can construct
            # an anonymous "formula" based on elements sorted by order number and
            # can do a = comparision to check if all elements are contained
            if len(quantities) > 1:
                raise Exception("HAS ONLY is not supported with zip")
            quantity = quantities[0]

            if quantity.has_only_quantity is None:
                raise Exception("HAS ONLY is not supported by %s" % quantity.name)

            def values():
                for predicates in predicate_zip_list:
                    if len(predicates) != 1:
                        raise Exception("Tuples not supported in HAS ONLY")
                    op, value = predicates[0]
                    if op != "=":
                        raise Exception("Predicated not supported in HAS ONLY")
                    if not isinstance(value, str):
                        raise Exception("Only strings supported in HAS ONLY")
                    yield value

            try:
                # Instead of the Optimade standard, the elements are combined
                # by the standard used by NOMAD.
                species, counts = get_hill_decomposition(list(values()))
                value = " ".join(["{}{}".format(s, "" if c == 1 else c) for s, c in zip(species, counts)])
            except KeyError:
                raise Exception("HAS ONLY is only supported for chemical symbols")

            return Q("term", **{quantity.has_only_quantity.name: value})
        else:
            raise NotImplementedError

        queries = [
            self._has_query(quantities, predicates) for predicates in predicate_zip_list
        ]
        args = {kind: queries}

        # Explicitly request that at least one term must match for the HAS ANY
        # query. ElasticSearch may either default to 0 or 1 depending on the
        # other search parameters:
        # https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-bool-query.html#bool-min-should-match
        if kind == "should":
            args["minimum_should_match"] = 1
        return Q("bool", **args)

    @v_args(inline=True)
    def property_first_comparison(self, quantity, query):
        # property_first_comparison: property *_rhs

        # If the quantity corresponds to a nested ES value, the query is
        # wrapped with the appropriate outer nested query.
        inner_query = query(quantity)
        if quantity.nested_path:
            return Q("nested", path=quantity.nested_path, query=inner_query)
        return inner_query

    @v_args(inline=True)
    def constant_first_comparison(self, value, op, quantity):
        # constant_first_comparison: constant OPERATOR ( non_string_value | ...not_implemented_string )
        if not isinstance(quantity, Quantity):
            raise Exception("Only quantities can be compared to constant values.")

        # If the quantity corresponds to a nested ES value, the query is
        # wrapped with the appropriate outer nested query.
        inner_query = self._query_op(quantity, _rev_cmp_operators[op], value)
        if quantity.nested_path:
            return Q("nested", path=quantity.nested_path, query=inner_query)
        return inner_query

    @v_args(inline=True)
    def value_op_rhs(self, op, value):
        # value_op_rhs: OPERATOR value
        return lambda quantity: self._query_op(
            quantity,
            op,
            value if quantity.converter is None else quantity.converter(value)
        )

    @v_args(inline=True)
    def boolean(self, value):
        if value == "TRUE":
            return True
        if value == "FALSE":
            return False
