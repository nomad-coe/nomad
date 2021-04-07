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

"""
API for retrieving material information.
"""
import re
import os
import math
import numpy as np
from typing import List
from collections import defaultdict

from flask_restplus import Resource, abort, fields, marshal
from flask import request, g
from elasticsearch_dsl import Search, Q, A, Text, Keyword, Boolean
from elasticsearch_dsl.query import Bool, Nested
from elasticsearch_dsl.utils import AttrDict
import ase.data
from lark import Lark

from nomad import config, infrastructure, search
from nomad.files import UploadFiles
from nomad.units import ureg
from nomad.atomutils import get_hill_decomposition
from nomad.datamodel.datamodel import EntryArchive
from nomad.datamodel.material import Material, Bulk, Method
from .materialtransformer import MElasticTransformer, MQuantity
from .api import api
from .auth import authenticate, create_authorization_predicate

ns = api.namespace("encyclopedia", description="Access materials data.")
missing_material_msg = "The specified material {} could not be retrieved. It either does not exists or requires authentication."
re_formula = re.compile(r"([A-Z][a-z]?)(\d*)")

# Parse using slightly modified Optimade grammar
with open("{}/encyclopedia_grammar.lark".format(os.path.dirname(__file__)), "r") as f:
    parser = Lark(f)


class MaterialAccessError(Exception):
    pass


class MaterialSearch():
    """Convenience class for material searches. Automatically ensures the
    correct visibility of materials when the search is constructed through the
    methods of his class.
    """
    def __init__(self):
        self._s = Search(index=config.elastic.materials_index_name)
        self._filters = []
        self._extra = {}
        self._q = None
        self.restricted = False

    def add_material_filter(self, query):
        """Adds material based filters.
        """
        self._filters.append(query)

    def add_material_aggregation(self, name, aggregation):
        """Adds material based aggregation.
        """
        self._s.aggs.bucket(name, aggregation)

    def add_calculation_filter(self, queries):
        """Adds calculation based filters. The visibility of calculations is
        automatically checked.
        """
        if not isinstance(queries, (list, tuple)):
            queries = [queries]
        nested_bool = Q(
            "bool",
            filter=queries,
        )
        nested_query = Q("nested", path="calculations", query=nested_bool)
        self._filters.append(nested_query)

    def includes(self, includes):
        self._extra["_source"] = {"includes": includes}

    def size(self, size):
        self._extra["size"] = size

    def extra(self, extra):
        self._extra = extra

    def s(self):
        # Wrap in bool query
        if self._q is None:
            query = Q("bool", filter=self._filters)
        else:
            query = self._q
        if not isinstance(query, Bool):
            query = Q("bool", filter=[query])

        # Split the query into material specific part and calculation specific
        # part. For now the queries are "separable", but this may not be the
        # case in the future...
        m_query = Q("bool")
        c_query = Q("bool")
        for q in query.must:
            c_query.must.append(q) if isinstance(q, Nested) else m_query.must.append(q)
        for q in query.filter:
            c_query.filter.append(q) if isinstance(q, Nested) else m_query.filter.append(q)
        for q in query.must_not:
            c_query.must_not.append(q) if isinstance(q, Nested) else m_query.must_not.append(q)
        for q in query.should:
            c_query.should.append(q) if isinstance(q, Nested) else m_query.should.append(q)

        # If restricted search is enabled, the order of nested/boolean queries
        # will be reversed depth-first in the calculation specific part.
        if self.restricted:
            def restrict(query):
                if isinstance(query, Bool):
                    # First restrict all inner queries, only after which the
                    # current query is restricted (=depth-first)
                    query.must = [q if isinstance(q, Nested) else restrict(q) for q in query.must]
                    query.filter = [q if isinstance(q, Nested) else restrict(q) for q in query.filter]
                    query.should = [q if isinstance(q, Nested) else restrict(q) for q in query.should]
                    query.must_not = [q if isinstance(q, Nested) else restrict(q) for q in query.must_not]

                    musts = [q.query if isinstance(q, Nested) else q for q in query.must]
                    filters = [q.query if isinstance(q, Nested) else q for q in query.filter]
                    shoulds = [q.query if isinstance(q, Nested) else q for q in query.should]
                    must_nots = [q.query if isinstance(q, Nested) else q for q in query.must_not]

                    inner_q = Q(
                        "bool",
                        filter=filters,
                        should=shoulds,
                        must=musts,
                        must_not=must_nots,
                    )
                    if len(shoulds) != 0:
                        inner_q.minimum_should_match = 1
                    outer_q = Q("nested", path="calculations", query=inner_q)
                    return outer_q
                else:
                    return query
            c_query = restrict(c_query)

        # Wrap in a boolean query if it is not already one.
        if not isinstance(c_query, Bool):
            c_query = Q("bool", filter=[c_query])

        # Merge calculation and material specific parts
        query = Q("bool")
        query.filter = c_query.filter + m_query.filter
        query.must = c_query.must + m_query.must
        query.must_not = c_query.must_not + m_query.must_not
        query.should = c_query.should + m_query.should

        # Add authentication filters on top of the query. This will make sure
        # that materials with only private calculations are excluded and that
        # private calculations are ignored in queries concerning calculations.
        def add_authentication(q):
            if isinstance(q, Nested):
                if q.path == "calculations":
                    nested_query = q.query
                    if not isinstance(nested_query, Bool):
                        nested_query = Q("bool", filter=[nested_query])
                    auth_filters = get_authentication_filters_material()
                    nested_query.filter.extend(auth_filters)
                    q.query = nested_query
            elif isinstance(q, Bool):
                for f in q.must:
                    add_authentication(f)
                for f in q.filter:
                    add_authentication(f)
                for f in q.should:
                    add_authentication(f)
                for f in q.must_not:
                    add_authentication(f)
        add_authentication(query)

        # This makes sure that materials with only private entries are always excluded
        query.filter.append(Q("nested", path="calculations", query=Q("bool", filter=get_authentication_filters_material())))

        # Enforce that all 'should' queries correspond to traditional 'or'
        # queries by enforcing minimum_should_match=1.
        def should_to_or(q):
            if isinstance(q, Bool):
                for f in q.must:
                    should_to_or(f)
                for f in q.filter:
                    should_to_or(f)
                for f in q.should:
                    should_to_or(f)
                for f in q.must_not:
                    should_to_or(f)
                if q.should:
                    q.minimum_should_match = 1
            if isinstance(q, Nested):
                should_to_or(q.query)
        should_to_or(query)

        s = self._s.query(query)
        # import json
        # print(json.dumps(s.to_dict(), indent=2))
        extra = self._extra
        s = s.extra(**extra)
        return s

    def execute(self):
        s = self.s()
        return s.execute()

    def calculations(self):
        """Executes the query and returns a list of visible calculations
        associated with the first found material. Currently fetches all
        calculations associated with a material. If the number of calculations
        per material increases significantly then the inner_hits available for
        nested queries should be used instead.

        Returns:
            List of visible calculations for the first material matching the
            constructed query.

        Raises:
            MaterialAccessError if the queried material could not be found.
        """
        source = self._extra.get("_source")
        if source is None:
            source = {}
            self._extra["_source"] = source
        includes = source.get("includes")
        if includes is None:
            includes = []
            source["includes"] = includes

        self._extra["_source"]["includes"].extend([
            "calculations.published",
            "calculations.with_embargo",
            "calculations.owners",
        ])
        response = self.execute()
        if response.hits.total == 0:
            raise MaterialAccessError

        material = response.hits[0]

        # Filter out calculations based on their visibility
        visible_calcs = []
        for calc in material.calculations:
            if calc.published and not calc.with_embargo:
                visible_calcs.append(calc)
            elif g.user is not None and g.user.user_id in calc.owners:
                visible_calcs.append(calc)
        return visible_calcs

    def from_query_string(self, query_string: str, restricted: bool) -> None:
        """Initializes this MaterialSearch instance from the specified search
        string.

        Args:
            query_string: The query string. E.g. 'crystal_system="cubic" AND
                material_type="bulk"'
            restricted: Whether the nested searches concerning calculations are
                combined or not.
        """
        try:
            parse_tree = parser.parse(query_string)
        except Exception as e:
            abort(400, message=(
                "Invalid query string: {}".format(e)
            ))

        # Transform parse tree into ES Query
        quantities: List[MQuantity] = [
            # Material level quantities
            MQuantity("elements", es_field="species", elastic_mapping_type=Text, has_only_quantity=MQuantity(name="species.keyword")),
            MQuantity("formula", es_field="species_and_counts", elastic_mapping_type=Text, has_only_quantity=MQuantity(name="species_and_counts.keyword")),
            MQuantity("material_id", es_field="material_id", elastic_mapping_type=Keyword),
            MQuantity("material_type", es_field="material_type", elastic_mapping_type=Keyword),
            MQuantity("material_name", es_field="material_name", elastic_mapping_type=Keyword),
            MQuantity("crystal_system", es_field="bulk.crystal_system", elastic_mapping_type=Keyword),
            MQuantity("space_group_number", es_field="bulk.space_group_number", elastic_mapping_type=Keyword),
            MQuantity("structure_type", es_field="bulk.structure_type", elastic_mapping_type=Keyword),
            # Calculation level (nested) quantities
            MQuantity("functional_type", es_field="calculations.method.functional_type", elastic_mapping_type=Keyword, nested_path="calculations"),
            MQuantity("basis_set", es_field="calculations.method.basis_set", elastic_mapping_type=Keyword, nested_path="calculations"),
            MQuantity("code_name", es_field="calculations.method.program_name", elastic_mapping_type=Keyword, nested_path="calculations"),
            MQuantity("has_band_structure", es_field="calculations.properties.has_electronic_band_structure", elastic_mapping_type=Boolean, nested_path="calculations"),
            MQuantity("has_dos", es_field="calculations.properties.has_electronic_dos", elastic_mapping_type=Boolean, nested_path="calculations"),
            MQuantity("has_thermal_properties", es_field="calculations.properties.has_thermodynamical_properties", elastic_mapping_type=Boolean, nested_path="calculations"),
            MQuantity("band_gap", es_field="calculations.properties.band_gap", elastic_mapping_type=Keyword, nested_path="calculations", converter=lambda x: (x * ureg.eV).to(ureg.joule).magnitude),
        ]
        transformer = MElasticTransformer(quantities)
        try:
            query = transformer.transform(parse_tree)
        except Exception as e:
            abort(400, message=(
                "Invalid query string: {}".format(e)
            ))

        self.restricted = restricted
        self._q = query


def get_authentication_filters_material():
    """Returns a list of queries that when placed inside a filter will leave
    out unpublished (of other users) or embargoed materials.
    """
    # Handle authentication
    filters = [Q('term', calculations__published=True), Q('term', calculations__with_embargo=False)]
    if g.user is not None and g.user.user_id is not None:
        q = filters[0] & filters[1]
        q = q | Q('term', calculations__owners=g.user.user_id)
        return [q]
    return filters


def get_authentication_filters_calc():
    """Returns a shared term filter that will leave out unpublished (of other
    users), embargoed or invalid entries in the calculations index.
    """
    # Handle authentication
    s = search.SearchRequest()
    if g.user is not None:
        s.owner('visible', user_id=g.user.user_id)
    else:
        s.owner('public')
    return [
        s.q,
        Q("term", encyclopedia__status="success"),
    ]


def get_range_filter(field, minimum=None, maximum=None, source_unit=None, target_unit=None):
    """For adding range filters
    """
    query_dict = {}
    if minimum is not None:
        if source_unit is None and target_unit is None:
            gte = minimum
        else:
            gte = (minimum * source_unit).to(target_unit).magnitude
        query_dict["gte"] = gte
    if maximum is not None:
        if source_unit is None and target_unit is None:
            lte = maximum
        else:
            lte = (maximum * source_unit).to(target_unit).magnitude
        query_dict["lte"] = lte
    query = Q("range", **{field: query_dict})
    return query


def rgetattr(obj, attr_name):
    """Used to perform attribute access based on a (possibly nested) attribute
    name given as string.
    """
    try:
        for attr in attr_name.split("."):
            obj = obj[attr]
    except KeyError:
        return None
    return obj


def get_es_doc_values(es_doc, mapping, keys=None):
    """Used to form a material definition for "materials/<material_id>" from
    the given ElasticSearch root document.
    """
    if keys is None:
        keys = mapping.keys()

    result = {}
    for key in keys:
        es_key = mapping[key]
        value = rgetattr(es_doc, es_key)
        if value is not None:
            result[key] = value

    return result


def read_archive(upload_id: str, calc_id: str) -> EntryArchive:
    """Used to read data from the archive.

    Args:
        upload_id: Upload id.
        calc_id: Calculation id.

    Returns:
        MSection: The section_run as MSection
        For each path, a dictionary containing the path as key and the returned
        section as value.
    """
    upload_files = UploadFiles.get(
        upload_id, is_authorized=create_authorization_predicate(upload_id, calc_id))

    with upload_files.read_archive(calc_id) as archive:
        data = archive[calc_id]
        root = EntryArchive.m_from_dict(data.to_dict())

    return root


def query_from_formula(formula: str) -> str:
    """Converts a formula into the corresponding element query.

    Args:
        formula: The requested formula.

    Returns:
        elasticsearch_dsl Query object.
    """
    element_list = []
    matches = re_formula.finditer(formula)

    prev_end = 0
    invalid = False
    for match in matches:
        if match.start() != prev_end:
            invalid = True
            break
        prev_end = match.end()
        groups = match.groups()
        symbol = groups[0]
        n_elem = groups[1]
        if symbol != "":
            if n_elem == "":
                element_list.append(symbol)
            else:
                element_list += [symbol] * int(n_elem)
    if prev_end != len(formula) or len(element_list) == 0:
        invalid = True
    if invalid:
        abort(400, message=(
            "Invalid chemical formula provided. Please use a single "
            "continuous string with capitalized element names and "
            "integer numbers for the element multiplicity. E.g. 'TiO2'"
        ))

    names, reduced_counts = get_hill_decomposition(element_list, reduced=True)
    query_list = []

    for name, count in zip(names, reduced_counts):
        if count == 1:
            query_list.append(name)
        else:
            query_list.append("{}{}".format(name, int(count)))
    query_string = " ".join(query_list)

    return query_string


def query_from_elements(elements: list) -> str:
    """Converts a list of elements into the corresponding element query.

    Args:
        elements: List of queried elmenents.

    Returns:
        elasticsearch_dsl Query object.
    """
    species, _ = get_hill_decomposition(elements)
    query_string = " ".join(species)

    return query_string


material_prop_map = {
    # General
    "material_id": "material_id",
    "formula": "formula",
    "formula_reduced": "formula_reduced",
    "material_type": "material_type",
    "material_name": "material_name",
    # Bulk
    "has_free_wyckoff_parameters": "bulk.has_free_wyckoff_parameters",
    "strukturbericht_designation": "bulk.strukturbericht_designation",
    "bravais_lattice": "bulk.bravais_lattice",
    "crystal_system": "bulk.crystal_system",
    "point_group": "bulk.point_group",
    "space_group_number": "bulk.space_group_number",
    "space_group_international_short_symbol": "bulk.space_group_international_short_symbol",
    "structure_type": "bulk.structure_type",
    "structure_prototype": "bulk.structure_prototype",
}
similarity = api.model("similarity", {
    # General
    "material_id": fields.String,
    "value": fields.Float,
    "formula": fields.String,
    "space_group_number": fields.Integer,
})
material_query = api.parser()
material_query.add_argument(
    "property",
    type=str,
    choices=tuple(material_prop_map.keys()),
    help="Optional single property to retrieve for the given material. If not specified, all properties will be returned.",
    location="args"
)
material_result = api.model("material_result", {
    # General
    "material_id": fields.String,
    "formula": fields.String,
    "formula_reduced": fields.String,
    "material_type": fields.String,
    "n_calculations": fields.Integer,
    # Bulk only
    "has_free_wyckoff_parameters": fields.Boolean,
    "strukturbericht_designation": fields.String,
    "material_name": fields.String,
    "bravais_lattice": fields.String,
    "crystal_system": fields.String,
    "point_group": fields.String,
    "space_group_number": fields.Integer,
    "space_group_international_short_symbol": fields.String,
    "structure_prototype": fields.String,
    "structure_type": fields.String,
    "similarity": fields.List(fields.Nested(similarity, skip_none=True), skip_none=True),
})


@ns.route("/materials/<string:material_id>")
class EncMaterialResource(Resource):
    @api.response(404, "The material does not exist")
    @api.response(200, "Metadata send", material_result)
    @api.expect(material_query)
    @api.marshal_with(material_result, skip_none=True)
    @api.doc("get_material", params={"material_id": "28 character identifier for the material."})
    @authenticate()
    def get(self, material_id):
        """Used to retrieve basic information related to a material.
        """
        # Parse request arguments
        args = material_query.parse_args()
        prop = args.get("property", None)
        if prop is not None:
            keys = [prop]
            es_keys = [material_prop_map[prop]]
        else:
            keys = list(material_prop_map.keys())
            es_keys = list(material_prop_map.values())

        # Get the material info, check that at least one calculation is visible
        s = MaterialSearch()
        s.add_material_filter(Q("term", material_id=material_id))
        s.includes(es_keys)
        response = s.execute()

        # No such material
        if response.hits.total == 0:
            abort(404, message=missing_material_msg.format(material_id))

        # Add values from ES entry
        entry = response[0]
        result = get_es_doc_values(entry, material_prop_map, keys)

        # Add similarity data that is stored in MongoDB.
        try:
            material = Material.m_def.a_mongo.get(material_id=material_id)
            dos_similarity = material.similarity.electronic_dos
        except KeyError:
            # No similarity data for this material
            pass
        else:
            # Only include similarity for materials that exist on the current
            # deployment to avoid dead links.
            similar_ids = dos_similarity.material_ids
            id_value_map = {key: value for key, value in zip(similar_ids, dos_similarity.values)}
            s = MaterialSearch()
            s.add_material_filter(Q("terms", material_id=similar_ids))
            s.includes(["material_id", "formula_reduced", "bulk.space_group_number"])
            s.size(5)
            response = s.execute()

            similarity = []
            for hit in response.hits:
                try:
                    similarity.append({
                        "material_id": hit.material_id,
                        "value": id_value_map[hit.material_id],
                        "formula": hit.formula_reduced,
                        "space_group_number": hit.bulk.space_group_number,
                    })
                except AttributeError:
                    pass
            if similarity:
                result["similarity"] = sorted(similarity, key=lambda x: x["value"], reverse=True)

        return result, 200


range_query = api.model("range_query", {
    "max": fields.Float(min=0),
    "min": fields.Float(min=0),
})
materials_query = api.model("materials_input", {
    "search_by": fields.Nested(api.model("search_query", {
        "exclusive": fields.Boolean(default=False, description="Set to True to enable exclusive element search."),
        "formula": fields.String(description="Chemical formula of the material as a string. Use a single continuous string with capitalized element names and integer numbers for the element multiplicity. The order of elements does not matter.", example="TiO2"),
        "elements": fields.List(fields.String(enum=ase.data.chemical_symbols[1:]), description="List of chemical species that the material should include. Use capitalized element name abbreviations.", example=["Ti", "O"]),
        "page": fields.Integer(default=1, min=1, description="Requested page number, indexing starts from 1.", example=1),
        "per_page": fields.Integer(default=25, min=1, description="Number of results per page.", example=10),
        "restricted": fields.Boolean(default=False, description="Select to restrict the query to individual calculations. If not selected, the query will combine results from several different calculations."),
    })),
    "query": fields.String(description="Single search string that supports combining any of the search parameters with AND, OR, NOT and parentheses."),
    "material_type": fields.List(fields.String(enum=list(Material.material_type.type)), description=Material.material_type.description),
    "material_name": fields.List(fields.String, description=Material.material_name.description),
    "structure_type": fields.List(fields.String, description=Bulk.structure_type.description),
    "space_group_number": fields.List(fields.Integer, min=1, max=230, description=Bulk.space_group_number.description),
    "crystal_system": fields.List(fields.String(enum=list(Bulk.crystal_system.type)), description=Bulk.crystal_system.description),
    "band_gap": fields.Nested(range_query, description="Band gap range in eV.", allow_null=True),
    "has_band_structure": fields.Boolean(description="Set to True if electronic band structure needs to be available for this material."),
    "has_dos": fields.Boolean(description="Set to True if electronic density of states needs to be available for this material."),
    "has_thermal_properties": fields.Boolean(description="Set to True if thermodynamical properties need to be available for this material."),
    "functional_type": fields.List(fields.String(enum=list(Method.functional_type.type)), description=Method.functional_type.description),
    "basis_set": fields.List(fields.String(enum=list(Method.basis_set.type)), description=Method.basis_set.description),
    "code_name": fields.List(fields.String(enum=list(Method.program_name.type)), description=Method.program_name.description),
})
pages_result = api.model("page_info", {
    "per_page": fields.Integer,
    "total": fields.Integer,
    "page": fields.Integer,
    "pages": fields.Integer,
})

materials_result = api.model("materials_result", {
    "total_results": fields.Integer(allow_null=False),
    "results": fields.List(fields.Nested(material_result, skip_none=True)),
    "pages": fields.Nested(pages_result, skip_none=True),
})


@ns.route("/materials/")
class EncMaterialsResource(Resource):
    @api.response(404, "No materials found")
    @api.response(400, "Bad request")
    @api.response(200, "OK", materials_result)
    @api.expect(materials_query, validate=False)
    @api.marshal_with(materials_result, skip_none=True)
    @api.doc("search_materials")
    @authenticate()
    def post(self):
        """Search materials based on their properties.
        """
        # Get query parameters as json
        try:
            data = marshal(request.get_json(), materials_query)
        except Exception as e:
            abort(400, message=str(e))

        s = MaterialSearch()
        query = data["query"]
        search_by = data["search_by"]
        exclusive = search_by["exclusive"]
        restricted = search_by["restricted"]

        # Initialize MaterialSearch from a query string
        if query is not None:
            s.from_query_string(query, restricted)
        # Initialize MaterialSearch from individual query terms
        else:
            # Material level filters
            if data["material_type"] is not None: s.add_material_filter(Q("terms", material_type=data["material_type"]))
            if data["material_name"] is not None: s.add_material_filter(Q("terms", material_name=data["material_name"]))
            if data["structure_type"] is not None: s.add_material_filter(Q("terms", bulk__structure_type=data["structure_type"]))
            if data["space_group_number"] is not None: s.add_material_filter(Q("terms", bulk__space_group_number=data["space_group_number"]))
            if data["crystal_system"] is not None: s.add_material_filter(Q("terms", bulk__crystal_system=data["crystal_system"]))

            # Calculation filters
            calc_filters = []
            if data["functional_type"] is not None: calc_filters.append(Q("terms", calculations__method__functional_type=data["functional_type"]))
            if data["basis_set"] is not None: calc_filters.append(Q("terms", calculations__method__basis_set=data["basis_set"]))
            if data["code_name"] is not None: calc_filters.append(Q("terms", calculations__method__program_name=data["code_name"]))
            if data["has_band_structure"] is not None: calc_filters.append(Q("term", calculations__properties__has_electronic_band_structure=data["has_band_structure"]))
            if data["has_dos"] is not None: calc_filters.append(Q("term", calculations__properties__has_electronic_dos=data["has_dos"]))
            if data["has_thermal_properties"] is not None: calc_filters.append(Q("term", calculations__properties__has_thermodynamical_properties=data["has_thermal_properties"]))
            if data["band_gap"] is not None: calc_filters.append(get_range_filter(
                "calculations.properties.band_gap",
                minimum=data["band_gap"].get("min"),
                maximum=data["band_gap"].get("max"),
                source_unit=ureg.eV,
                target_unit=ureg.J,
            ))
            if restricted:
                s.add_calculation_filter(calc_filters)
            else:
                for f in calc_filters:
                    s.add_calculation_filter(f)

            # The given list of species/formula is reformatted with the Hill system into a
            # query string. With exclusive search we look for exact match, with
            # non-exclusive search we look for match that includes at least all
            # species, possibly even more.
            formula = search_by["formula"]
            elements = search_by["elements"]
            if formula is not None or elements is not None:
                if formula is not None:
                    query_string = query_from_formula(formula)
                    path = "species_and_counts"
                elif elements is not None:
                    query_string = query_from_elements(elements)
                    path = "species"
                if exclusive:
                    s.add_material_filter(Q("term", **{"{}.keyword".format(path): query_string}))
                else:
                    s.add_material_filter(Q(
                        "match",
                        **{path: {"query": query_string, "operator": "and"}},
                    ))

        # Execute query
        page = search_by["page"]
        per_page = search_by["per_page"]
        s.extra({
            "size": per_page,
            "from": (page - 1) * per_page,
            "sort": [{"formula_reduced": {"order": "asc"}}],
            "_source": {"includes": list(material_prop_map.values())},
        })
        response = s.execute()

        # Form final response
        pages = {
            "page": page,
            "per_page": per_page,
            "pages": math.ceil(response.hits.total / per_page),
            "total": response.hits.total,
        }

        # Gather the number of visible calculation for each returned material
        # with an aggregation
        if len(response) != 0:
            material_ids = [x.material_id for x in response]
            s2 = MaterialSearch()
            s2.size(0)
            matched = s2._s.aggs.bucket("matched", A("filter", filter=Q("terms", material_id=material_ids)))
            materials = matched.bucket("materials", A("terms", field="material_id", size=len(material_ids)))
            nested = materials.bucket("nested", A("nested", path="calculations"))
            nested.bucket(
                "visible",
                A("filter", filter=Q("bool", filter=get_authentication_filters_material()))
            )
            response2 = s2.execute()
            agg_dict = {}
            for agg in response2.aggs.matched.materials:
                agg_dict[agg.key] = agg.nested.visible.doc_count

        # Form the final list of results
        result_list = []
        for x in response:
            res = get_es_doc_values(x, material_prop_map, list(material_prop_map.keys()))
            material_id = x.material_id
            res["n_calculations"] = agg_dict[material_id]
            result_list.append(res)

        return {"results": result_list, "pages": pages}, 200


groups_result = api.model("groups_result", {
    "groups_eos": fields.Raw,
    "groups_par": fields.Raw,
})


@ns.route("/materials/<string:material_id>/groups")
class EncGroupsResource(Resource):
    @api.response(404, "Material not found")
    @api.response(400, "Bad request")
    @api.response(200, "OK", groups_result)
    @api.marshal_with(groups_result)
    @api.doc("get_material_groups", params={"material_id": "28 character identifier for the material."})
    @authenticate()
    def get(self, material_id):
        """Returns a summary of the calculation groups that were identified for this material.

        Two types of groups are reported: equation of state groups and
        parameter variation groups. Equation of state groups contain
        calculations with identical method and material, but different volume.
        Parameter variation groups contain identical structure but different
        methods. The response contains dictionaries for both groups
        ('groups_eos' and 'groups_par'). These dictionaries map a group id with
        a list of calculation ids.
        """
        # Get full entry for this material
        s = MaterialSearch()
        s.add_material_filter(Q("term", material_id=material_id))
        s.extra({
            "_source": {"includes": [
                "calculations.calc_id",
                "calculations.method.group_eos_id",
                "calculations.method.group_parametervariation_id",
                "calculations.properties.energies.energy_total",
                "calculations.idealized_structure.cell_volume",
            ]},
            "size": 1,
        })

        # Raise error if material not found
        try:
            calculations = s.calculations()
        except MaterialAccessError:
            abort(404, message=missing_material_msg.format(material_id))

        groups_eos = defaultdict(list)
        groups_param = defaultdict(list)
        for calc in calculations:
            try:
                calc.properties.energies.energy_total
                calc.idealized_structure.cell_volume
            except AttributeError:
                continue
            try:
                group_eos_id = calc.method.group_eos_id
                if group_eos_id:
                    groups_eos[group_eos_id].append(calc.calc_id)
            except AttributeError:
                pass
            try:
                group_param_id = calc.method.group_parametervariation_id
                if group_param_id:
                    groups_param[group_param_id].append(calc.calc_id)
            except AttributeError:
                pass

        # Filter out groups with too few entries
        for key, items in list(groups_eos.items()):
            if len(items) < 4:
                del groups_eos[key]
        for key, items in list(groups_param.items()):
            if len(items) < 2:
                del groups_param[key]

        # Return results
        result = {
            "groups_eos": groups_eos,
            "groups_par": groups_param,
        }

        return result, 200


group_result = api.model("group_result", {
    "calculations": fields.List(fields.String, description="List of calculation ids."),
    "energies": fields.List(fields.Float, description="List of total energies."),
    "volumes": fields.List(fields.Float, description="List of cell volumes."),
})


@ns.route("/materials/<string:material_id>/groups/<string:group_type>/<string:group_id>")
class EncGroupResource(Resource):
    @api.response(404, "Group not found")
    @api.response(400, "Bad request")
    @api.response(200, "OK", group_result)
    @api.marshal_with(group_result)
    @api.doc("get_material_group", params={
        "material_id": "28 character identifier for the material.",
        "group_type": "Type of group. Valid options are: 'eos' and 'par'.",
        "group_id": "28 character identifier for the group.",
    })
    @authenticate()
    def get(self, material_id, group_type, group_id):
        """Used to query detailed information about a specific calculation group.
        """
        # Find entries for the given material, which have EOS or parameter
        # variation hashes set.
        if group_type == "eos":
            group_id_source = "group_eos_id"
        elif group_type == "par":
            group_id_source = "group_parametervariation_id"
        else:
            abort(400, message="Unsupported group type.")

        s = MaterialSearch()
        s.add_material_filter(Q("term", material_id=material_id))
        s.extra({
            "_source": {"includes": [
                "calculations.calc_id",
                "calculations.properties.energies.energy_total",
                "calculations.idealized_structure.cell_volume",
                "calculations.method." + group_id_source,
            ]},
            "size": 1,
        })

        # Raise error if material not found
        try:
            calculations = s.calculations()
        except MaterialAccessError:
            abort(404, message=missing_material_msg.format(material_id))

        # Gather groups from the calculations
        calcs = []
        energies = []
        volumes = []
        for calc in calculations:
            try:
                i_group_id = getattr(calc.method, group_id_source)
                if i_group_id == group_id:
                    calcs.append(calc.calc_id)
                    volumes.append(calc.idealized_structure.cell_volume)
                    energies.append(calc.properties.energies.energy_total)
            except Exception:
                pass

        # Sort results by energy
        energies = np.array(energies)
        volumes = np.array(volumes)
        calcs = np.array(calcs)
        order = energies.argsort()
        energies = energies[order]
        volumes = volumes[order]
        calcs = calcs[order]

        # Return results
        group_dict = {
            "calculations": calcs.tolist(),
            "energies": energies.tolist(),
            "volumes": volumes.tolist(),
        }

        return group_dict, 200


calc_prop_map = {
    "calc_id": "calc_id",
    "upload_id": "upload_id",
    "code_name": "method.program_name",
    "code_version": "method.program_version",
    "functional_type": "method.functional_type",
    "basis_set_type": "method.basis_set",
    "core_electron_treatment": "method.core_electron_treatment",
    "run_type": "workflow.workflow_type",
    "has_dos": "properties.has_electronic_dos",
    "has_band_structure": "properties.has_electronic_band_structure",
    "has_thermal_properties": "properties.has_thermodynamical_properties",
}
calculation_result = api.model("calculation_result", {
    "calc_id": fields.String,
    "upload_id": fields.String,
    "code_name": fields.String,
    "code_version": fields.String,
    "functional_type": fields.String,
    "basis_set_type": fields.String,
    "core_electron_treatment": fields.String(default="unavailable"),
    "run_type": fields.String(default="unavailable"),
    "has_dos": fields.Boolean,
    "has_band_structure": fields.Boolean,
    "has_thermal_properties": fields.Boolean,
})
representatives_result = api.model("representatives_result", {
    "idealized_structure": fields.String,
    "electronic_band_structure": fields.String,
    "electronic_dos": fields.String,
    "thermodynamical_properties": fields.String,
})
calculations_result = api.model("calculations_result", {
    "total_results": fields.Integer,
    "results": fields.List(fields.Nested(calculation_result)),
    "representatives": fields.Nested(representatives_result, skip_none=True),
})


@ns.route("/materials/<string:material_id>/calculations")
class EncCalculationsResource(Resource):
    @api.response(404, "Material not found")
    @api.response(400, "Bad request")
    @api.response(200, "OK", calculations_result)
    @api.doc("get_material_calculations")
    @api.marshal_with(calculations_result)
    @authenticate()
    def get(self, material_id):
        """Used to return information about all calculations related to the given material.

        Returns a list of all calculations and a representative calculation for
        few select quantities that are shown in the material overview page.
        """
        s = MaterialSearch()
        s.add_material_filter(Q("term", material_id=material_id))
        s.extra({"_source": {"includes": ["calculations"]}})

        def calc_score(entry):
            """Custom scoring function used to sort results by their
            "quality". Currently built to mimic the scoring that was used
            in the old Encyclopedia GUI. Primarily sorts by quality measure,
            ties are broken by alphabetic sorting of entry_id in order to
            return consistent results.
            """
            score = 0
            functional_score = {
                "GGA": 100
            }
            code_score = {
                "VASP": 3,  # Prefer VASP data as it is the "cleanest" on average
                "FHI-aims": 2,
                "Quantum Espresso": 1,
            }
            code_name = entry.method.program_name
            functional = entry.method.functional_type
            try:
                has_bs = entry.properties.has_electronic_band_structure
            except AttributeError:
                has_bs = False
            try:
                has_dos = entry.properties.has_electronic_dos
            except AttributeError:
                has_dos = False
            score += functional_score.get(functional, 0)
            score += code_score.get(code_name, 0)
            if has_dos and has_bs:
                score += 10

            return (score, entry.calc_id)

        # Raise error if material not found
        try:
            calculations = s.calculations()
        except MaterialAccessError:
            abort(404, message=missing_material_msg.format(material_id))

        # Sort calculations by "quality"
        sorted_calc = sorted(calculations, key=lambda x: calc_score(x), reverse=True)

        # Get the requested representative properties
        representatives = {}
        representatives["idealized_structure"] = sorted_calc[0].calc_id
        thermo_found = False
        bs_found = False
        dos_found = False
        for calc in sorted_calc:
            if not hasattr(calc, "properties"):
                continue

            if not thermo_found and calc.properties.has_thermodynamical_properties:
                representatives["thermodynamical_properties"] = calc.calc_id
                thermo_found = True
            if not bs_found and calc.properties.has_electronic_band_structure:
                representatives["electronic_band_structure"] = calc.calc_id
                bs_found = True
            if not dos_found and calc.properties.has_electronic_dos:
                representatives["electronic_dos"] = calc.calc_id
                dos_found = True
            if thermo_found and bs_found and dos_found:
                break

        # Create result JSON
        results = []
        for entry in sorted_calc:
            calc_dict = get_es_doc_values(entry, calc_prop_map)
            results.append(calc_dict)

        result = {
            "total_results": len(results),
            "results": results,
            "representatives": representatives,
        }

        return result, 200


histogram = api.model("histogram", {
    "occurrences": fields.List(fields.Integer),
    "values": fields.List(fields.Float),
})
statistics_query = api.model("statistics_query", {
    "calculations": fields.List(fields.String),
    "properties": fields.List(fields.String),
    "n_histogram_bins": fields.Integer,
})
statistics = api.model("statistics", {
    "min": fields.Float,
    "max": fields.Float,
    "avg": fields.Float,
    "histogram": fields.Nested(histogram, skip_none=True)
})
statistics_result = api.model("statistics_result", {
    "cell_volume": fields.Nested(statistics, skip_none=True),
    "atomic_density": fields.Nested(statistics, skip_none=True),
    "mass_density": fields.Nested(statistics, skip_none=True),
    "lattice_a": fields.Nested(statistics, skip_none=True),
    "lattice_b": fields.Nested(statistics, skip_none=True),
    "lattice_c": fields.Nested(statistics, skip_none=True),
    "alpha": fields.Nested(statistics, skip_none=True),
    "beta": fields.Nested(statistics, skip_none=True),
    "gamma": fields.Nested(statistics, skip_none=True),
    "band_gap": fields.Nested(statistics, skip_none=True),
})
property_map = {
    "cell_volume": "encyclopedia.material.idealized_structure.cell_volume",
    "atomic_density": "encyclopedia.properties.atomic_density",
    "mass_density": "encyclopedia.properties.mass_density",
    "lattice_a": "encyclopedia.material.idealized_structure.lattice_parameters.a",
    "lattice_b": "encyclopedia.material.idealized_structure.lattice_parameters.b",
    "lattice_c": "encyclopedia.material.idealized_structure.lattice_parameters.c",
    "alpha": "encyclopedia.material.idealized_structure.lattice_parameters.alpha",
    "beta": "encyclopedia.material.idealized_structure.lattice_parameters.beta",
    "gamma": "encyclopedia.material.idealized_structure.lattice_parameters.gamma",
    "band_gap": "encyclopedia.properties.band_gap",
}


@ns.route("/materials/<string:material_id>/statistics")
class EncStatisticsResource(Resource):
    @api.response(404, "Suggestion not found")
    @api.response(400, "Bad request")
    @api.response(200, "OK", statistics_result)
    @api.expect(statistics_query, validate=False)
    @api.marshal_with(statistics_result, skip_none=True)
    @api.doc("get_material_statistics", params={"material_id": "28 character identifier for the material."})
    @authenticate()
    def post(self, material_id):
        """Used to return statistics related to the specified material and
        calculations.
        """
        # Get query parameters as json
        try:
            data = marshal(request.get_json(), statistics_query)
        except Exception as e:
            abort(400, message=str(e))

        # Find entries for the given material.
        bool_query = Q(
            "bool",
            filter=get_authentication_filters_calc() + [
                Q("term", encyclopedia__material__material_id=material_id),
                Q("terms", calc_id=data["calculations"]),
            ]
        )

        s = Search(index=config.elastic.index_name)
        s = s.query(bool_query)
        s = s.extra(**{
            "size": 0,
        })

        # Add statistics aggregations for each requested property
        properties = data["properties"]
        for prop in properties:
            stats_agg = A("stats", field=property_map[prop])
            s.aggs.bucket("{}_stats".format(prop), stats_agg)

        # No hits on the top query level
        response = s.execute()
        if response.hits.total == 0:
            abort(404, message="The given calculations could not be found for material {}".format(material_id))

        # Run a second query that creates histograms with fixed size buckets
        # based on the min and max from previous query. Might make more sense
        # to use the mean and sigma to define the range?
        s = Search(index=config.elastic.index_name)
        s = s.query(bool_query)
        s = s.extra(**{
            "size": 0,
        })
        n_bins = data["n_histogram_bins"]
        for prop in properties:
            stats = getattr(response.aggs, "{}_stats".format(prop))
            if stats.count == 0:
                continue
            interval = (stats.max * 1.001 - stats.min) / n_bins
            if interval == 0:
                interval = 1
            hist_agg = A("histogram", field=property_map[prop], interval=interval, offset=stats.min, min_doc_count=0)
            s.aggs.bucket("{}_hist".format(prop), hist_agg)
        response_hist = s.execute()

        # Return results
        result = {}
        for prop in properties:
            stats = getattr(response.aggs, "{}_stats".format(prop))
            if stats.count == 0:
                continue
            hist = getattr(response_hist.aggs, "{}_hist".format(prop))
            occurrences = [x.doc_count for x in hist.buckets]
            values = [x.key for x in hist.buckets]
            result[prop] = {
                "min": stats.min,
                "max": stats.max,
                "avg": stats.avg,
                "histogram": {
                    "occurrences": occurrences,
                    "values": values,
                }
            }

        return result, 200


wyckoff_variables_result = api.model("wyckoff_variables_result", {
    "x": fields.Float,
    "y": fields.Float,
    "z": fields.Float,
})
wyckoff_set_result = api.model("wyckoff_set_result", {
    "wyckoff_letter": fields.String,
    "indices": fields.List(fields.Integer),
    "element": fields.String,
    "variables": fields.Nested(wyckoff_variables_result, skip_none=True),
})
lattice_parameters = api.model("lattice_parameters", {
    "a": fields.Float,
    "b": fields.Float,
    "c": fields.Float,
    "alpha": fields.Float,
    "beta": fields.Float,
    "gamma": fields.Float,
})

idealized_structure_result = api.model("idealized_structure_result", {
    "atom_labels": fields.List(fields.String),
    "atom_positions": fields.List(fields.List(fields.Float)),
    "lattice_vectors": fields.List(fields.List(fields.Float)),
    "lattice_vectors_primitive": fields.List(fields.List(fields.Float)),
    "lattice_parameters": fields.Nested(lattice_parameters, skip_none=True),
    "periodicity": fields.List(fields.Boolean),
    "number_of_atoms": fields.Integer,
    "cell_volume": fields.Float,
    "wyckoff_sets": fields.List(fields.Nested(wyckoff_set_result, skip_none=True)),
})

calculation_property_map = {
    "lattice_parameters": {
        "source": "es",
        "path": "encyclopedia.material.idealized_structure.lattice_parameters"
    },
    "energies": {
        "source": "es",
        "path": "encyclopedia.properties.energies",
    },
    "mass_density": {
        "source": "es",
        "path": "encyclopedia.properties.mass_density",
    },
    "atomic_density": {
        "source": "es",
        "path": "encyclopedia.properties.atomic_density",
    },
    "cell_volume": {
        "source": "es",
        "path": "encyclopedia.material.idealized_structure.cell_volume"
    },
    "band_gap": {
        "source": "es",
        "path": "encyclopedia.properties.band_gap"
    },
    "electronic_band_structure": {
        "source": "es",
        "path": "encyclopedia.properties.electronic_band_structure"
    },
    "electronic_dos": {
        "source": "es",
        "path": "encyclopedia.properties.electronic_dos"
    },
    "phonon_band_structure": {
        "source": "es",
        "path": "encyclopedia.properties.phonon_band_structure"
    },
    "phonon_dos": {
        "source": "es",
        "path": "encyclopedia.properties.phonon_dos"
    },
    "thermodynamical_properties": {
        "source": "es",
        "path": "encyclopedia.properties.thermodynamical_properties"
    },
    "wyckoff_sets": {
        "source": "archive",
        "path": "section_metadata/encyclopedia/material/idealized_structure/wyckoff_sets"
    },
    "idealized_structure": {
        "source": "archive",
        "path": "section_metadata/encyclopedia/material/idealized_structure"
    },
}

calculation_property_query = api.model("calculation_query", {
    "properties": fields.List(fields.String(enum=list(calculation_property_map.keys())), description="List of calculation properties to return."),
})
energies = api.model("energies", {
    "energy_total": fields.Float,
    "energy_total_T0": fields.Float,
    "energy_free": fields.Float,
})
electronic_band_structure = api.model("electronic_band_structure", {
    "reciprocal_cell": fields.List(fields.List(fields.Float)),
    "brillouin_zone": fields.Raw,
    "section_k_band_segment": fields.Raw,
    "section_band_gap": fields.Raw,
})
electronic_dos = api.model("electronic_dos", {
    "dos_energies": fields.List(fields.Float),
    "dos_values": fields.List(fields.List(fields.Float)),
})
calculation_property_result = api.model("calculation_property_result", {
    "lattice_parameters": fields.Nested(lattice_parameters, skip_none=True),
    "energies": fields.Nested(energies, skip_none=True),
    "mass_density": fields.Float,
    "atomic_density": fields.Float,
    "cell_volume": fields.Float,
    "wyckoff_sets": fields.Nested(wyckoff_set_result, skip_none=True),
    "idealized_structure": fields.Nested(idealized_structure_result, skip_none=True),
    "band_gap": fields.Float,
    "electronic_band_structure": fields.Nested(electronic_band_structure, skip_none=True),
    "electronic_dos": fields.Nested(electronic_dos, skip_none=True),
    "phonon_band_structure": fields.Raw,
    "phonon_dos": fields.Raw,
    "thermodynamical_properties": fields.Raw,
})


@ns.route("/materials/<string:material_id>/calculations/<string:calc_id>")
class EncCalculationResource(Resource):
    @api.response(404, "Material or calculation not found")
    @api.response(400, "Bad request")
    @api.response(200, "OK", calculation_property_result)
    @api.expect(calculation_property_query, validate=False)
    @api.marshal_with(calculation_property_result, skip_none=True)
    @api.doc("get_calculation")
    @authenticate()
    def post(self, material_id, calc_id):
        """Get properties from a specific calculation related to a material.
        """
        # Get query parameters as json
        try:
            data = marshal(request.get_json(), calculation_property_query)
        except Exception as e:
            abort(400, message=str(e))

        s = Search(index=config.elastic.index_name)
        query = Q(
            "bool",
            filter=get_authentication_filters_calc() + [
                Q("term", encyclopedia__material__material_id=material_id),
                Q("term", calc_id=calc_id),
            ]
        )
        s = s.query(query)

        # Create dictionaries for requested properties
        references = []
        properties = data["properties"]
        es_properties = {}
        mongo_properties = {}
        arch_properties = {}
        ref_properties = set((
            "electronic_dos",
            "electronic_band_structure",
            "thermodynamical_properties",
            "phonon_dos",
            "phonon_band_structure",
        ))
        for prop in properties:
            source = calculation_property_map[prop]["source"]
            path = calculation_property_map[prop]["path"]
            if source == "es":
                es_properties[prop] = path
                if prop in ref_properties:
                    references.append(prop)
            elif source == "mongo":
                mongo_properties[prop] = path
            elif source == "archive":
                arch_properties[prop] = path

        # The query is filtered already on the ES side so we don't need to
        # transfer so much data.
        sources = [
            "upload_id",
            "calc_id",
            "encyclopedia",
        ]
        sources += list(es_properties.values())

        s = s.extra(**{
            "_source": {"includes": sources},
            "size": 1,
        })

        response = s.execute()

        # No such material
        if len(response) == 0:
            abort(404, message=(
                "Could not retrieve calculation {} for material {}. The "
                "entry either does not exist or requires authentication."
                .format(calc_id, material_id))
            )

        # Add references that are to be read from the archive
        for ref in references:
            arch_path = response[0]
            arch_path = rgetattr(arch_path, es_properties[ref])
            if arch_path is not None:
                arch_properties[ref] = arch_path
            del es_properties[ref]

        # If any of the requested properties require data from the Archive, the
        # file is opened and read.
        result = {}
        if len(arch_properties) != 0:
            entry = response[0]
            upload_id = entry.upload_id
            calc_id = entry.calc_id
            root = read_archive(
                upload_id,
                calc_id,
            )

            # Add results from archive
            for key, arch_path in arch_properties.items():
                try:
                    value = root[arch_path]

                    # Replace unnormalized thermodynamical properties with
                    # normalized ones and turn into dict
                    if key == "thermodynamical_properties":
                        specific_heat_capacity = value.specific_heat_capacity.magnitude.tolist()
                        specific_free_energy = value.specific_vibrational_free_energy_at_constant_volume.magnitude.tolist()
                        specific_heat_capacity = [x if np.isfinite(x) else None for x in specific_heat_capacity]
                        specific_free_energy = [x if np.isfinite(x) else None for x in specific_free_energy]
                    if isinstance(value, list):
                        value = [x.m_to_dict() for x in value]
                    else:
                        value = value.m_to_dict()
                    if key == "thermodynamical_properties":
                        del value["thermodynamical_property_heat_capacity_C_v"]
                        del value["vibrational_free_energy_at_constant_volume"]
                        value["specific_heat_capacity"] = specific_heat_capacity
                        value["specific_vibrational_free_energy_at_constant_volume"] = specific_free_energy

                    # DOS results are simplified.
                    if key == "electronic_dos":
                        if "dos_energies_normalized" in value:
                            value["dos_energies"] = value["dos_energies_normalized"]
                            del value["dos_energies_normalized"]
                        if "dos_values_normalized" in value:
                            value["dos_values"] = value["dos_values_normalized"]
                            del value["dos_values_normalized"]

                    # Pre-calculate k-path length to be used as x-coordinate in
                    # plots. If the VBM and CBM information is needed later, it
                    # can be added as indices along the path. The exact k-points
                    # and occupations are removed to save some bandwidth.
                    if key == "electronic_band_structure" or key == "phonon_band_structure":
                        segments = value["section_k_band_segment"]
                        k_path_length = 0
                        for segment in segments:
                            k_points = np.array(segment["band_k_points"])
                            segment_length = np.linalg.norm(k_points[-1, :] - k_points[0, :])
                            k_path_distances = k_path_length + np.linalg.norm(k_points - k_points[0, :], axis=1)
                            k_path_length += segment_length
                            segment["k_path_distances"] = k_path_distances.tolist()
                            del segment["band_k_points"]
                            if "band_occupations" in segment:
                                del segment["band_occupations"]

                    result[key] = value
                except (AttributeError, KeyError):
                    abort(500, "Could not find the requested resource.")

        # Add results from ES
        for prop, es_source in es_properties.items():
            value = rgetattr(response[0], es_source)
            if value is not None:
                if isinstance(value, AttrDict):
                    value = value.to_dict()
                result[prop] = value

        # Add results from Mongo
        if len(mongo_properties) != 0:
            mongo_db = infrastructure.mongo_client[config.mongo.db_name]
            archives = mongo_db['archive']
            archive = archives.find_one({"_id": calc_id})
            for prop, mongo_source in mongo_properties.items():
                value = rgetattr(archive, mongo_source)
                if value is not None:
                    result[prop] = value

        return result, 200


suggestions_map = {
    "code_name": "dft.code_name",
    "basis_set": "dft.basis_set",
    "functional_type": "encyclopedia.method.functional_type",
    "structure_type": "bulk.structure_type",
    "strukturbericht_designation": "bulk.strukturbericht_designation",
}
suggestions_query = api.parser()
suggestions_query.add_argument(
    "property",
    type=str,
    choices=("code_name", "structure_type", "strukturbericht_designation", "basis_set", "functional_type"),
    help="The property name for which suggestions are returned.",
    location="args"
)
suggestions_result = api.model("suggestions_result", {
    "code_name": fields.List(fields.String),
    "basis_set": fields.List(fields.String),
    "functional_type": fields.List(fields.String),
    "structure_type": fields.List(fields.String),
    "strukturbericht_designation": fields.List(fields.String),
})


@ns.route("/suggestions")
class EncSuggestionsResource(Resource):
    @api.response(404, "Suggestion not found")
    @api.response(400, "Bad request")
    @api.response(200, "OK", suggestions_result)
    @api.expect(suggestions_query, validate=False)
    @api.marshal_with(suggestions_result, skip_none=True)
    @api.doc("get_material_suggestions")
    @authenticate()
    def get(self):
        """Dynamically retrieves a list of unique values for the given property.
        """
        # Uses terms aggregation to return all unique terms for the requested
        # field. Without using composite aggregations there is a size limit for
        # the number of aggregation buckets. This should, however, not be a
        # problem since the number of unique values is low for all supported
        # properties.

        # Parse request arguments
        args = suggestions_query.parse_args()
        prop = args.get("property", None)

        # Material level suggestions
        if prop in {"structure_type", "strukturbericht_designation"}:
            s = MaterialSearch()
            s.size(0)
            s.add_material_aggregation("suggestions", A("terms", field=suggestions_map[prop], size=999))
        # Calculation level suggestions
        elif prop in {"code_name", "basis_set", "functional_type"}:
            s = Search(index=config.elastic.index_name)
            query = Q(
                "bool",
                filter=get_authentication_filters_calc()
            )
            s = s.query(query)
            s = s.extra(**{
                "size": 0,
            })

            terms_agg = A("terms", field=suggestions_map[prop], size=999)
            s.aggs.bucket("suggestions", terms_agg)
        else:
            raise ValueError("No suggestion available for '{}'".format(prop))

        # Gather unique values into a list
        response = s.execute()
        suggestions = [x.key for x in response.aggs.suggestions.buckets]

        return {prop: suggestions}, 200


report_query = api.model("report_query", {
    "server": fields.String,
    "username": fields.String,
    "email": fields.String,
    "first_name": fields.String,
    "last_name": fields.String,
    "category": fields.String,
    "subcategory": fields.String(allow_null=True),
    "representatives": fields.Raw(Raw=True),
    "message": fields.String,
})


@ns.route("/materials/<string:material_id>/reports")
class ReportsResource(Resource):
    @api.response(500, "Error sending report")
    @api.response(400, "Bad request")
    @api.response(204, "Report succesfully sent")
    @api.expect(report_query)
    @api.doc("post_material_report", params={"material_id": "28 character identifier for the material."})
    @authenticate(required=True)
    def post(self, material_id):
        """Post an error report on a material. Requires authentication.
        """
        # Get query parameters as json
        try:
            query = marshal(request.get_json(), report_query)
        except Exception as e:
            abort(400, message=str(e))

        # Send the report as an email
        query["material_id"] = material_id
        representatives = query["representatives"]
        if representatives is not None:
            representatives = "\n" + "\n".join(["  {}: {}".format(key, value) for key, value in representatives.items()])
            query["representatives"] = representatives
        mail = (
            "Server: {server}\n\n"
            "Username: {username}\n"
            "First name: {first_name}\n"
            "Last name: {last_name}\n"
            "Email: {email}\n\n"
            "Material id: {material_id}\n"
            "Category: {category}\n"
            "Subcategory: {subcategory}\n"
            "Representative calculations: {representatives}\n\n"
            "Message: {message}"
        ).format(**query)
        try:
            infrastructure.send_mail(
                name="webmaster", email="support@nomad-lab.eu", message=mail, subject='Encyclopedia error report')
        except Exception as e:
            abort(500, message="Error sending error report email.")
        return "", 204
