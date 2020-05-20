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

"""
The encyclopedia API of the nomad@FAIRDI APIs.
"""
import re
import math

from flask_restplus import Resource, abort, fields, marshal
from flask import request
from elasticsearch_dsl import Search, Q, A

from nomad import config
from nomad.units import ureg
from nomad.atomutils import get_hill_decomposition
from .api import api

ns = api.namespace("encyclopedia", description="Access encyclopedia metadata.")
re_formula = re.compile(r"([A-Z][a-z]?)(\d*)")

material_prop_map = {
    # General
    "material_id": "encyclopedia.material.material_id",
    "formula": "encyclopedia.material.formula",
    "formula_reduced": "encyclopedia.material.formula_reduced",
    "system_type": "encyclopedia.material.material_type",
    # Bulk only
    "has_free_wyckoff_parameters": "encyclopedia.material.bulk.has_free_wyckoff_parameters",
    "strukturbericht_designation": "encyclopedia.material.bulk.strukturbericht_designation",
    "material_name": "encyclopedia.material.material_name",
    "bravais_lattice": "encyclopedia.material.bulk.bravais_lattice",
    "crystal_system": "encyclopedia.material.bulk.crystal_system",
    "point_group": "encyclopedia.material.bulk.point_group",
    "space_group_number": "encyclopedia.material.bulk.space_group_number",
    "space_group_international_short_symbol": "encyclopedia.material.bulk.space_group_international_short_symbol",
    "structure_prototype": "encyclopedia.material.bulk.structure_prototype",
    "structure_type": "encyclopedia.material.bulk.structure_type",
}


def get_material(es_doc, keys):
    """Used to form a material definition for "materials/<material_id>" from
    the given ElasticSearch root document.
    """
    result = {}
    for key in keys:
        es_key = material_prop_map[key]
        try:
            value = es_doc
            for part in es_key.split("."):
                value = getattr(value, part)
        except AttributeError:
            value = None
        result[key] = value

    return result


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
    "system_type": fields.String,
    # Bulk only
    "has_free_wyckoff_parameters": fields.String,
    "strukturbericht_designation": fields.String,
    "material_name": fields.String,
    "bravais_lattice": fields.String,
    "crystal_system": fields.String,
    "point_group": fields.String,
    "space_group_number": fields.Integer,
    "space_group_international_short_symbol": fields.String,
    "structure_prototype": fields.String,
    "structure_type": fields.String,
})


@ns.route("/materials/<string:material_id>")
class EncMaterialResource(Resource):
    @api.response(404, "The material does not exist")
    @api.response(200, "Metadata send", fields.Raw)
    @api.doc("material/<material_id>")
    @api.expect(material_query)
    @api.marshal_with(material_result, skip_none=True)
    def get(self, material_id):
        """Used to retrive basic information related to the specified material.
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

        # Find the first public entry with this material id and take
        # information from there. In principle all other entries should have
        # the same information.
        s = Search(index=config.elastic.index_name)

        # Since we are looking for an exact match, we use filter context
        # together with term search for speed (instead of query context and
        # match search)
        query = Q(
            "bool",
            filter=[
                Q("term", published=True),
                Q("term", with_embargo=False),
                Q("term", encyclopedia__material__material_id=material_id),
            ]
        )
        s = s.query(query)

        # The query is collapsed already on the ES side so we don"t need to
        # transfer so much data.
        s = s.extra(**{
            "collapse": {"field": "encyclopedia.material.material_id"},
            "_source": {"includes": es_keys},
        })

        response = s.execute()

        # No such material
        if len(response) == 0:
            abort(404, message="There is no material {}".format(material_id))

        # Create result JSON
        entry = response[0]
        result = get_material(entry, keys)

        return result, 200


range_query = api.model("range_query", {
    "max": fields.Float,
    "min": fields.Float,
})
materials_query = api.model("materials_input", {
    "search_by": fields.Nested(api.model("search_query", {
        "exclusive": fields.Boolean(default=False),
        "formula": fields.String,
        "element": fields.String,
        "page": fields.Integer(default=1),
        "per_page": fields.Integer(default=25),
        "pagination": fields.Boolean,
        "mode": fields.String(default="collapse"),
    })),
    "material_name": fields.List(fields.String),
    "structure_type": fields.List(fields.String),
    "space_group_number": fields.List(fields.Integer),
    "system_type": fields.List(fields.String),
    "crystal_system": fields.List(fields.String),
    "band_gap": fields.Nested(range_query, description="Band gap range in eV."),
    "band_gap_direct": fields.Boolean,
    "has_band_structure": fields.Boolean,
    "has_dos": fields.Boolean,
    "has_fermi_surface": fields.Boolean,
    "has_thermal_properties": fields.Boolean,
    "functional_type": fields.List(fields.String),
    "basis_set_type": fields.List(fields.String),
    "code_name": fields.List(fields.String),
    "mass_density": fields.Nested(range_query, description="Mass density range in kg / m ** 3."),
})
materials_result = api.model("materials_result", {
    "total_results": fields.Integer(allow_null=False),
    "results": fields.List(fields.Nested(material_result)),
    "pages": fields.Nested(api.model("page_info", {
        "per_page": fields.Integer,
        "total": fields.Integer,
        "page": fields.Integer,
        "pages": fields.Integer,
    })),
    "es_query": fields.String(allow_null=False),
})


@ns.route("/materials")
class EncMaterialsResource(Resource):
    @api.response(404, "No materials found")
    @api.response(400, "Bad request")
    @api.response(200, "Metadata send", fields.Raw)
    @api.expect(materials_query, validate=False)
    @api.marshal_with(materials_result, skip_none=True)
    @api.doc("materials")
    def post(self):
        """Used to query a list of materials with the given search options.
        """
        # Get query parameters as json
        try:
            data = marshal(request.get_json(), materials_query)
        except Exception as e:
            abort(400, message=str(e))

        filters = []
        must_nots = []
        musts = []

        # Add term filters
        filters.append(Q("term", published=True))
        filters.append(Q("term", with_embargo=False))

        def add_terms_filter(source, target, query_type="terms"):
            if data[source]:
                filters.append(Q(query_type, **{target: data[source]}))

        add_terms_filter("material_name", "encyclopedia.material.material_name")
        add_terms_filter("structure_type", "encyclopedia.material.bulk.structure_type")
        add_terms_filter("space_group_number", "encyclopedia.material.bulk.space_group_number")
        add_terms_filter("system_type", "encyclopedia.material.material_type")
        add_terms_filter("crystal_system", "encyclopedia.material.bulk.crystal_system")
        add_terms_filter("band_gap_direct", "encyclopedia.properties.band_gap_direct", query_type="term")
        add_terms_filter("functional_type", "encyclopedia.method.functional_type")
        add_terms_filter("basis_set_type", "dft.basis_set")
        add_terms_filter("code_name", "dft.code_name")

        # Add exists filters
        def add_exists_filter(source, target):
            param = data[source]
            if param is not None:
                query = Q("exists", field=target)
                if param is True:
                    filters.append(query)
                elif param is False:
                    must_nots.append(query)

        add_exists_filter("has_thermal_properties", "encyclopedia.properties.thermodynamical_properties")
        add_exists_filter("has_band_structure", "encyclopedia.properties.electronic_band_structure")
        add_exists_filter("has_dos", "encyclopedia.properties.electronic_dos")
        add_exists_filter("has_fermi_surface", "encyclopedia.properties.fermi_surface")

        # Add range filters
        def add_range_filter(source, target, source_unit=None, target_unit=None):
            param = data[source]
            query_dict = {}
            if param["min"] is not None:
                if source_unit is None and target_unit is None:
                    gte = param["min"]
                else:
                    gte = (param["min"] * source_unit).to(target_unit).magnitude
                query_dict["gte"] = gte
            if param["max"] is not None:
                if source_unit is None and target_unit is None:
                    lte = param["max"]
                else:
                    lte = (param["max"] * source_unit).to(target_unit).magnitude
                query_dict["lte"] = lte
            if len(query_dict) != 0:
                query = Q("range", **{target: query_dict})
                filters.append(query)

        add_range_filter("band_gap", "encyclopedia.properties.band_gap", ureg.eV, ureg.J)
        add_range_filter("mass_density", "encyclopedia.properties.mass_density")

        # Create query for elements or formula
        search_by = data["search_by"]
        mode = search_by["mode"]
        formula = search_by["formula"]
        elements = search_by["element"]
        exclusive = search_by["exclusive"]

        if formula is not None:
            # Here we determine a list of atom types. The types may occur
            # multiple times and at multiple places.
            element_list = []
            matches = re_formula.finditer(formula)
            for match in matches:
                groups = match.groups()
                symbol = groups[0]
                count = groups[1]
                if symbol != "":
                    if count == "":
                        element_list.append(symbol)
                    else:
                        element_list += [symbol] * int(count)

            # The given list of species is reformatted with the Hill system
            # into a query string. The counts are reduced by the greatest
            # common divisor.
            names, reduced_counts = get_hill_decomposition(element_list, reduced=True)
            query_string = []
            for name, count in zip(names, reduced_counts):
                if count == 1:
                    query_string.append(name)
                else:
                    query_string.append("{}{}".format(name, int(count)))
            query_string = " ".join(query_string)

            # With exclusive search we look for exact match
            if exclusive:
                filters.append(Q("term", **{"encyclopedia.material.species_and_counts.keyword": query_string}))
            # With non-exclusive search we look for match that includes at
            # least all parts of the formula, possibly even more.
            else:
                musts.append(Q(
                    "match",
                    encyclopedia__material__species_and_counts={"query": query_string, "operator": "and"}
                ))
        elif elements is not None:
            # The given list of species is reformatted with the Hill system into a query string
            species, _ = get_hill_decomposition(elements.split(","))
            query_string = " ".join(species)

            # With exclusive search we look for exact match
            if exclusive:
                filters.append(Q("term", **{"encyclopedia.material.species.keyword": query_string}))
            # With non-exclusive search we look for match that includes at
            # least all species, possibly even more.
            else:
                musts.append(Q(
                    "match",
                    encyclopedia__material__species={"query": query_string, "operator": "and"}
                ))

        page = search_by["page"]
        per_page = search_by["per_page"]
        bool_query = Q(
            "bool",
            filter=filters,
            must_not=must_nots,
            must=musts,
        )

        # 1: The paginated approach: No way to know the amount of matches,
        # but can return aggregation results in a quick fashion including
        # the number of matches entries per material.
        if mode == "aggregate":
            after = None
            # The loop is awkward, but emulates the old behaviour until the GUI is adapted.
            for _ in range(page):

                # The top query filters out entries based on the user query
                s = Search(index=config.elastic.index_name)
                s = s.query(bool_query)

                # The materials are grouped by using three aggregations:
                # "Composite" to enable scrolling, "Terms" to enable selecting
                # by material_id and "Top Hits" to fetch a single
                # representative material document. Unnecessary fields are
                # filtered to reduce data transfer.
                terms_agg = A("terms", field="encyclopedia.material.material_id")
                composite_kwargs = {"sources": {"materials": terms_agg}, "size": per_page}
                if after is not None:
                    composite_kwargs["after"] = after
                composite_agg = A("composite", **composite_kwargs)
                composite_agg.metric("representative", A(
                    "top_hits",
                    size=1,
                    _source={"includes": list(material_prop_map.values())},
                ))
                s.aggs.bucket("materials", composite_agg)

                # We ignore the top level hits
                s = s.extra(**{
                    "size": 0,
                })

                response = s.execute()
                materials = response.aggs.materials.buckets
                if len(materials) == 0:
                    abort(404, message="No materials found for the given search criteria or pagination.")
                after = response.aggs.materials["after_key"]

            # Gather results from aggregations
            result_list = []
            materials = response.aggs.materials.buckets
            keys = list(material_prop_map.keys())
            for material in materials:
                representative = material["representative"][0]
                mat_dict = get_material(representative, keys)
                mat_dict["n_of_calculations"] = material.doc_count
                result_list.append(mat_dict)

            # Page information is incomplete for aggregations
            pages = {
                "page": page,
                "per_page": per_page,
            }
        # 2. Collapse approach. Quickly provides a list of materials
        # corresponding to the query, offers full pagination, doesn"t include
        # the number of matches per material.
        elif mode == "collapse":
            s = Search(index=config.elastic.index_name)
            s = s.query(bool_query)
            s = s.extra(**{
                "collapse": {"field": "encyclopedia.material.material_id"},
                "size": per_page,
                "from": (page - 1) * per_page,
            })

            # Execute query
            response = s.execute()

            # No matches
            if len(response) == 0:
                abort(404, message="No materials found for the given search criteria or pagination.")

            # Loop over materials
            result_list = []
            keys = list(material_prop_map.keys())
            for material in response:
                mat_result = get_material(material, keys)
                result_list.append(mat_result)

            # Full page information available for collapse
            pages = {
                "page": page,
                "per_page": per_page,
                "pages": math.ceil(response.hits.total / per_page),
                "total": response.hits.total,
            }

        result = {
            "results": result_list,
            "total_results": len(result_list),
            "es_query": s.to_dict(),
            "pages": pages,
        }
        return result, 200


group_result = api.model("group_result", {
    "calculation_list": fields.List(fields.String),
    "energy_minimum": fields.Float,
    "group_hash": fields.String,
    "group_type": fields.String,
    "nr_of_calculations": fields.Integer,
    "representative_calculation_id": fields.String,
})
groups_result = api.model("groups_result", {
    "total_groups": fields.Integer(allow_null=False),
    "groups": fields.List(fields.Nested(group_result)),
})
group_source = {
    "includes": [
        "calc_id",
        "encyclopedia.properties.energies.energy_total",
    ]
}


@ns.route("/materials/<string:material_id>/groups")
class EncGroupsResource(Resource):
    @api.response(404, "Material not found")
    @api.response(400, "Bad request")
    @api.response(200, "Metadata send", fields.Raw)
    @api.expect(material_query, validate=False)
    @api.marshal_with(groups_result, skip_none=True)
    @api.doc("enc_materials")
    def get(self, material_id):

        # Find entries for the given material, which have EOS or parameter
        # variation hashes set.
        bool_query = Q(
            "bool",
            filter=[
                Q("term", published=True),
                Q("term", with_embargo=False),
                Q("term", encyclopedia__material__material_id=material_id),
            ],
            must=[
                Q("exists", field="encyclopedia.properties.energies.energy_total"),
            ],
            should=[
                Q("exists", field="encyclopedia.method.group_eos_hash"),
                Q("exists", field="encyclopedia.method.group_parametervariation_hash"),
            ],
            minimum_should_match=1,  # At least one of the should query must match
        )

        s = Search(index=config.elastic.index_name)
        s = s.query(bool_query)

        # Bucket the calculations by the group hashes. Only create a bucket if an
        # above-minimum number of documents are found.
        group_eos_bucket = A("terms", field="encyclopedia.method.group_eos_hash", min_doc_count=4)
        group_param_bucket = A("terms", field="encyclopedia.method.group_parametervariation_hash", min_doc_count=2)

        # calc_id and energy should be extracted for each matched document. The
        # documents are sorted by energy so that the minimum energy one can be
        # easily extracted. A maximum request size is set in order to limit the
        # result size. ES also has an index-level property
        # "index.max_inner_result_window" that limits the number of results
        # that an inner result can contain.
        energy_aggregation = A(
            "top_hits",
            _source=group_source,
            sort=[{"encyclopedia.properties.energies.energy_total": {"order": "asc"}}],
            size=100,
        )
        group_eos_bucket.bucket("energies", energy_aggregation)
        group_param_bucket.bucket("energies", energy_aggregation)
        s.aggs.bucket("groups_eos", group_eos_bucket)
        s.aggs.bucket("groups_param", group_param_bucket)

        # We ignore the top level hits
        s = s.extra(**{
            "size": 0,
        })

        # No hits on the top query level
        response = s.execute()
        n_hits = response.hits.total
        if n_hits == 0:
            abort(404, message="The specified material could not be found.")

        # Collect information for each group from the aggregations
        groups = []
        groups_eos = response.aggs.groups_eos.buckets
        groups_param = response.aggs.groups_param.buckets

        def get_group(group, group_type, group_hash):
            hits = group.energies.hits
            calculations = [doc.calc_id for doc in hits]
            group_dict = {
                "group_hash": group_hash,
                "group_type": group_type,
                "nr_of_calculations": len(calculations),
                "representative_calculation_id": hits[0].calc_id,
                "calculation_list": calculations,
                "energy_minimum": hits[0].encyclopedia.properties.energies.energy_total,
            }
            return group_dict

        for group in groups_eos:
            groups.append(get_group(group, "equation of state", group.key))
        for group in groups_param:
            groups.append(get_group(group, "parameter variation", group.key))

        # Return results
        result = {
            "groups": groups,
            "total_groups": len(groups),
        }
        return result, 200
