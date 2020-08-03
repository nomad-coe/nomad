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
import numpy as np

from flask_restplus import Resource, abort, fields, marshal
from flask import request
from elasticsearch_dsl import Search, Q, A
from elasticsearch_dsl.utils import AttrDict

from nomad import config, files, infrastructure
from nomad.units import ureg
from nomad.atomutils import get_hill_decomposition
from nomad.datamodel.datamodel import EntryArchive
from .api import api
from .auth import authenticate

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
        result[key] = value

    return result


def get_enc_filter():
    """Returns a shared term filter that will leave out unpublished, embargoed
    or invalid entries.
    """
    return [
        Q("term", published=True),
        Q("term", with_embargo=False),
        Q("term", encyclopedia__status="success"),
    ]


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
    "n_matches": fields.Integer,
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
})


@ns.route("/materials/<string:material_id>")
class EncMaterialResource(Resource):
    @api.response(404, "The material does not exist")
    @api.response(200, "Metadata send", fields.Raw)
    @api.doc("material/<material_id>")
    @api.expect(material_query)
    @api.marshal_with(material_result, skip_none=True)
    def get(self, material_id):
        """Used to retrieve basic information related to the specified
        material.
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
        query = Q(
            "bool",
            filter=get_enc_filter() + [
                Q("term", encyclopedia__material__material_id=material_id),
            ]
        )
        s = s.query(query)

        # Only one representative entry is returned by collapsing the results.
        s = s.extra(**{
            "_source": {"includes": es_keys},
            "size": 1,
            "collapse": {"field": "encyclopedia.material.material_id"},
        })
        response = s.execute()

        # No such material
        if len(response) == 0:
            abort(404, message="There is no material {}".format(material_id))

        # Add values from ES entry
        entry = response[0]
        result = get_es_doc_values(entry, material_prop_map, keys)

        return result, 200


range_query = api.model("range_query", {
    "max": fields.Float,
    "min": fields.Float,
})
materials_after = api.model("materials_after", {
    "materials": fields.String,
})
materials_query = api.model("materials_input", {
    "search_by": fields.Nested(api.model("search_query", {
        "exclusive": fields.Boolean(default=False),
        "formula": fields.String,
        "element": fields.String,
        "page": fields.Integer(default=1),
        "after": fields.Nested(materials_after, allow_null=True),
        "per_page": fields.Integer(default=25),
        "pagination": fields.Boolean,
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
pages_result = api.model("page_info", {
    "per_page": fields.Integer,
    "total": fields.Integer,
    "page": fields.Integer,
    "pages": fields.Integer,
    "after": fields.Nested(materials_after),
})

materials_result = api.model("materials_result", {
    "total_results": fields.Integer(allow_null=False),
    "results": fields.List(fields.Nested(material_result, skip_none=True)),
    "pages": fields.Nested(pages_result, skip_none=True),
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

        # The queries that correspond to AND queries typically need to access
        # multiple calculations at once to find the material ids that
        # correspond to the query. To implement this behaviour we need to run
        # an initial aggregation that checks that the requested properties are
        # present for a material. This is a a very crude solution that does not
        # scale to complex queries, but I'm not sure we can do much better
        # until we have a separate index for materials.
        property_map = {
            "has_thermal_properties": "encyclopedia.properties.thermodynamical_properties",
            "has_band_structure": "encyclopedia.properties.electronic_band_structure",
            "has_dos": "encyclopedia.properties.electronic_dos",
            "has_fermi_surface": "encyclopedia.properties.fermi_surface",
        }
        requested_properties = []
        # The size is set very large because all the results need to be
        # returned. We cannot get the results in a paginated way with composite
        # aggregation, because pipeline aggregations are not compatible with
        # them.
        agg_parent = A("terms", field="encyclopedia.material.material_id", size=5000000)
        for key, value in property_map.items():
            if data[key] is True:
                agg = A("filter", exists={"field": value})
                agg_parent.bucket(key, agg)
                requested_properties.append(key)
        if len(requested_properties) > 1:
            bool_query = Q(
                "bool",
                filter=get_enc_filter(),
            )
            s = Search(index=config.elastic.index_name)
            s = s.query(bool_query)
            s.aggs.bucket("materials", agg_parent)
            buckets_path = {x: "{}._count".format(x) for x in requested_properties}
            script = " && ".join(["params.{} > 0".format(x) for x in requested_properties])
            agg_parent.pipeline("selector", A(
                "bucket_selector",
                buckets_path=buckets_path,
                script=script,
            ))
            s = s.extra(**{
                "size": 0,
            })
            response = s.execute()
            material_ids = [x["key"] for x in response.aggs.materials.buckets]
            if len(material_ids) == 0:
                abort(404, message="No materials found for the given search criteria or pagination.")

        # After finding the material ids that fill the AND conditions, continue
        # with a simple OR query.
        filters = get_enc_filter()
        must_nots = []
        musts = []

        def add_terms_filter(source, target, query_type="terms"):
            if data[source] is not None:
                filters.append(Q(query_type, **{target: data[source]}))

        if len(requested_properties) > 1:
            filters.append(Q("terms", encyclopedia__material__material_id=material_ids))
        add_terms_filter("material_name", "encyclopedia.material.material_name")
        add_terms_filter("structure_type", "encyclopedia.material.bulk.structure_type")
        add_terms_filter("space_group_number", "encyclopedia.material.bulk.space_group_number")
        add_terms_filter("system_type", "encyclopedia.material.material_type")
        add_terms_filter("crystal_system", "encyclopedia.material.bulk.crystal_system")
        add_terms_filter("band_gap_direct", "encyclopedia.properties.band_gap_direct", query_type="term")
        add_terms_filter("functional_type", "encyclopedia.method.functional_type")
        add_terms_filter("basis_set_type", "dft.basis_set")
        add_terms_filter("code_name", "dft.code_name")

        # Add exists filters if only one property was requested. The initial
        # aggregation will handlei multiple simultaneous properties.
        def add_exists_filter(source, target):
            param = data[source]
            if param is not None:
                query = Q("exists", field=target)
                if param is True:
                    filters.append(query)
                elif param is False:
                    must_nots.append(query)
        if len(requested_properties) == 1:
            prop_name = requested_properties[0]
            add_exists_filter(prop_name, property_map[prop_name])

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
        after = search_by["after"]
        bool_query = Q(
            "bool",
            filter=filters,
            must_not=must_nots,
            must=musts,
        )

        # The top query filters out entries based on the user query
        s = Search(index=config.elastic.index_name)
        s = s.query(bool_query)

        # 1: The paginated approach: No way to know the amount of materials,
        # but can return aggregation results in a quick fashion including
        # the number of calculation entries per material.
        mode = "collapse"
        if mode == "aggregation":
            # The materials are grouped by using three aggregations:
            # "Composite" to enable scrolling, "Terms" to enable selecting
            # by material_id and "Top Hits" to fetch a single
            # representative material document. Unnecessary fields are
            # filtered to reduce data transfer.
            terms_agg = A("terms", field="encyclopedia.material.material_id")
            composite_kwargs = {"sources": {"materials": terms_agg}, "size": per_page}

            # The number of matched materials is only requested on the first
            # search, not for each page.
            if after is not None:
                composite_kwargs["after"] = after
            else:
                cardinality_agg = A("cardinality", field="encyclopedia.material.material_id", precision_threshold=1000)
                s.aggs.metric("n_materials", cardinality_agg)

            composite_agg = A("composite", **composite_kwargs)
            composite_agg.metric("representative", A(
                "top_hits",
                size=1,
                _source={"includes": list(material_prop_map.values())},
            ))
            s.aggs.bucket("materials", composite_agg)

            # We ignore the top level hits and sort by reduced material formula.
            s = s.extra(**{
                "size": 0,
            })

            response = s.execute()
            materials = response.aggs.materials.buckets
            if len(materials) == 0:
                abort(404, message="No materials found for the given search criteria or pagination.")
            after_new = response.aggs.materials["after_key"]

            # Gather results from aggregations
            result_list = []
            materials = response.aggs.materials.buckets
            keys = list(material_prop_map.keys())
            for material in materials:
                representative = material["representative"][0]
                mat_dict = get_es_doc_values(representative, material_prop_map, keys)
                mat_dict["n_matches"] = material.doc_count
                result_list.append(mat_dict)

            # Page information is incomplete for aggregations
            pages = {
                "page": page,
                "per_page": per_page,
                "after": after_new,
            }
            if after is None:
                n_materials = response.aggs.n_materials.value
                pages["total"] = n_materials

        # 2. Collapse approach. Quickly provides a list of materials
        # corresponding to the query, offers full pagination, the number of
        # matches per material needs to be requested with a separate query.
        elif mode == "collapse":
            s = Search(index=config.elastic.index_name)
            s = s.query(bool_query)

            # Add cardinality aggregation that gives out the total number of materials
            cardinality_agg = A("cardinality", field="encyclopedia.material.material_id", precision_threshold=1000)
            s.aggs.metric("n_materials", cardinality_agg)

            s = s.extra(**{
                "collapse": {"field": "encyclopedia.material.material_id"},
                "size": per_page,
                "from": (page - 1) * per_page,
                "sort": [{"encyclopedia.material.formula_reduced": {"order": "asc"}}],
                "explain": True,
            })

            # Execute query
            response = s.execute()

            # No matches
            if len(response) == 0:
                abort(404, message="No materials found for the given search criteria or pagination.")

            # Gather number of entries per material with a separate query
            material_ids = [x.encyclopedia.material.material_id for x in response]
            s = Search(index=config.elastic.index_name)
            bool_query = Q(
                "bool",
                filter=Q("terms", encyclopedia__material__material_id=material_ids),
            )
            s2 = s.query(bool_query)
            s2.aggs.bucket("n_matches", A("terms", field="encyclopedia.material.material_id"))
            response2 = s2.execute()
            matmap = {x.key: x.doc_count for x in response2.aggs.n_matches}

            # Loop over materials
            result_list = []
            keys = list(material_prop_map.keys())
            for material in response:
                # Get values from the collapsed doc
                mat_result = get_es_doc_values(material, material_prop_map, keys)
                mat_id = material.encyclopedia.material.material_id
                mat_result["n_matches"] = matmap[mat_id]
                result_list.append(mat_result)

            # Full page information available for collapse
            pages = {
                "page": page,
                "per_page": per_page,
                "pages": math.ceil(response.hits.total / per_page),
                "total": response.aggs.n_materials.value,
            }

        result = {
            "results": result_list,
            "pages": pages,
        }
        return result, 200


groups_result = api.model("groups_result", {
    "groups_eos": fields.Raw,
    "groups_par": fields.Raw,
})


@ns.route("/materials/<string:material_id>/groups")
class EncGroupsResource(Resource):
    @api.response(404, "Material not found")
    @api.response(400, "Bad request")
    @api.response(200, "Metadata send", fields.Raw)
    @api.marshal_with(groups_result)
    @api.doc("enc_materials")
    def get(self, material_id):
        """Returns a summary of the calculation groups that were identified for
        this material.
        """
        # Find entries for the given material, which have EOS or parameter
        # variation hashes set.
        bool_query = Q(
            "bool",
            filter=get_enc_filter() + [Q("term", encyclopedia__material__material_id=material_id)],
            must=[
                Q("exists", field="encyclopedia.properties.energies.energy_total"),
                Q("exists", field="encyclopedia.material.idealized_structure.cell_volume"),
            ],
            should=[
                Q("exists", field="encyclopedia.method.group_eos_id"),
                Q("exists", field="encyclopedia.method.group_parametervariation_id"),
            ],
            minimum_should_match=1,  # At least one of the should query must match
        )

        s = Search(index=config.elastic.index_name)
        s = s.query(bool_query)

        # Bucket the calculations by the group hashes. Only create a bucket if an
        # above-minimum number of documents are found.
        group_eos_bucket = A("terms", field="encyclopedia.method.group_eos_id", min_doc_count=4)
        group_param_bucket = A("terms", field="encyclopedia.method.group_parametervariation_id", min_doc_count=2)
        calc_aggregation = A(
            "top_hits",
            _source={"includes": ["calc_id"]},
            sort=[{"encyclopedia.properties.energies.energy_total": {"order": "asc"}}],
            size=100,
        )
        group_eos_bucket.bucket("calculations", calc_aggregation)
        group_param_bucket.bucket("calculations", calc_aggregation)
        s.aggs.bucket("groups_eos", group_eos_bucket)
        s.aggs.bucket("groups_param", group_param_bucket)

        # We ignore the top level hits
        s = s.extra(**{
            "size": 0,
        })

        # Collect information for each group from the aggregations
        response = s.execute()
        groups_eos = {group.key: [calc.calc_id for calc in group.calculations.hits] for group in response.aggs.groups_eos.buckets}
        groups_param = {group.key: [calc.calc_id for calc in group.calculations.hits] for group in response.aggs.groups_param.buckets}

        # Return results
        result = {
            "groups_eos": groups_eos,
            "groups_par": groups_param,
        }

        return result, 200


group_result = api.model("group_result", {
    "calculations": fields.List(fields.String),
    "energies": fields.List(fields.Float),
    "volumes": fields.List(fields.Float),
})
group_source = {
    "includes": [
        "calc_id",
        "encyclopedia.properties.energies.energy_total",
        "encyclopedia.material.idealized_structure.cell_volume",
    ]
}


@ns.route("/materials/<string:material_id>/groups/<string:group_type>/<string:group_id>")
class EncGroupResource(Resource):
    @api.response(404, "Group not found")
    @api.response(400, "Bad request")
    @api.response(200, "Metadata send", fields.Raw)
    @api.marshal_with(group_result)
    @api.doc("enc_group")
    def get(self, material_id, group_type, group_id):
        """Used to query detailed information for a specific calculation group.
        """
        # Find entries for the given material, which have EOS or parameter
        # variation hashes set.
        if group_type == "eos":
            group_id_source = "encyclopedia.method.group_eos_id"
        elif group_type == "par":
            group_id_source = "encyclopedia.method.group_parametervariation_id"
        else:
            abort(400, message="Unsupported group type.")

        bool_query = Q(
            "bool",
            filter=get_enc_filter() + [
                Q("term", encyclopedia__material__material_id=material_id),
                Q("term", **{group_id_source: group_id}),
            ],
        )

        s = Search(index=config.elastic.index_name)
        s = s.query(bool_query)

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
        s.aggs.bucket("groups_eos", energy_aggregation)

        # We ignore the top level hits
        s = s.extra(**{
            "size": 0,
        })

        # Collect information for each group from the aggregations
        response = s.execute()

        hits = response.aggs.groups_eos.hits
        calculations = [doc.calc_id for doc in hits]
        energies = [doc.encyclopedia.properties.energies.energy_total for doc in hits]
        volumes = [doc.encyclopedia.material.idealized_structure.cell_volume for doc in hits]
        group_dict = {
            "calculations": calculations,
            "energies": energies,
            "volumes": volumes,
        }

        return group_dict, 200


suggestions_map = {
    "code_name": "dft.code_name",
    "structure_type": "encyclopedia.material.bulk.structure_type",
}
suggestions_query = api.parser()
suggestions_query.add_argument(
    "property",
    type=str,
    choices=("code_name", "structure_type"),
    help="The property name for which suggestions are returned.",
    location="args"
)
suggestions_result = api.model("suggestions_result", {
    "code_name": fields.List(fields.String),
    "structure_type": fields.List(fields.String),
})


@ns.route("/suggestions")
class EncSuggestionsResource(Resource):
    @api.response(404, "Suggestion not found")
    @api.response(400, "Bad request")
    @api.response(200, "Metadata send", fields.Raw)
    @api.expect(suggestions_query, validate=False)
    @api.marshal_with(suggestions_result, skip_none=True)
    @api.doc("enc_suggestions")
    def get(self):

        # Parse request arguments
        args = suggestions_query.parse_args()
        prop = args.get("property", None)

        # Use aggregation to return all unique terms for the requested field.
        # Without using composite aggregations there is a size limit for the
        # number of aggregation buckets. This should, however, not be a problem
        # since the number of unique values is low for all supported properties.
        s = Search(index=config.elastic.index_name)
        query = Q(
            "bool",
            filter=get_enc_filter()
        )
        s = s.query(query)
        s = s.extra(**{
            "size": 0,
        })

        terms_agg = A("terms", field=suggestions_map[prop])
        s.aggs.bucket("suggestions", terms_agg)

        # Gather unique values into a list
        response = s.execute()
        suggestions = [x.key for x in response.aggs.suggestions.buckets]

        return {prop: suggestions}, 200


calc_prop_map = {
    "calc_id": "calc_id",
    "code_name": "dft.code_name",
    "code_version": "dft.code_version",
    "functional_type": "encyclopedia.method.functional_type",
    "basis_set_type": "dft.basis_set",
    "core_electron_treatment": "encyclopedia.method.core_electron_treatment",
    "run_type": "encyclopedia.calculation.calculation_type",
    "has_dos": "encyclopedia.properties.electronic_dos",
    "has_band_structure": "encyclopedia.properties.electronic_band_structure",
    "has_thermal_properties": "encyclopedia.properties.thermodynamical_properties",
    "has_phonon_dos": "encyclopedia.properties.phonon_dos",
    "has_phonon_band_structure": "encyclopedia.properties.phonon_band_structure",
}
calculation_result = api.model("calculation_result", {
    "calc_id": fields.String,
    "code_name": fields.String,
    "code_version": fields.String,
    "functional_type": fields.String,
    "basis_set_type": fields.String,
    "core_electron_treatment": fields.String,
    "run_type": fields.String,
    "has_dos": fields.Boolean,
    "has_band_structure": fields.Boolean,
    "has_thermal_properties": fields.Boolean,
    "has_phonon_dos": fields.Boolean,
    "has_phonon_band_structure": fields.Boolean,
})
representatives_result = api.model("representatives_result", {
    "idealized_structure": fields.String,
    "electronic_band_structure": fields.String,
    "electronic_dos": fields.String,
    "thermodynamical_properties": fields.String,
})
calculations_result = api.model("calculations_result", {
    "total_results": fields.Integer,
    "pages": fields.Nested(pages_result),
    "results": fields.List(fields.Nested(calculation_result)),
    "representatives": fields.Nested(representatives_result, skip_none=True),
})


@ns.route("/materials/<string:material_id>/calculations")
class EncCalculationsResource(Resource):
    @api.response(404, "Suggestion not found")
    @api.response(400, "Bad request")
    @api.response(200, "Metadata send", fields.Raw)
    @api.doc("enc_calculations")
    def get(self, material_id):
        """Used to return all calculations related to the given material. Also
        returns a representative calculation for each property shown in the
        overview page.
        """
        s = Search(index=config.elastic.index_name)
        query = Q(
            "bool",
            filter=get_enc_filter() + [
                Q("term", encyclopedia__material__material_id=material_id),
            ]
        )
        s = s.query(query)

        # The query is filtered already on the ES side so we don"t need to
        # transfer so much data.
        s = s.extra(**{
            "_source": {"includes": list(calc_prop_map.values()) + ["dft.xc_functional"]},
            "size": 10000,
            "from": 0,
        })
        response = s.execute()

        # No such material
        if len(response) == 0:
            abort(404, message="There is no material {}".format(material_id))

        # Add representative properties. It might be possible to write a custom
        # ES scoring mechanism or aggregation to also perform the selection.
        representatives = {}

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
                "FHI-aims": 3,
                "VASP": 2,
                "Quantum Espresso": 1,
            }
            code_name = entry.dft.code_name
            functional = entry.dft.xc_functional
            has_dos = rgetattr(entry, "encyclopedia.properties.electronic_band_structure") is not None
            has_bs = rgetattr(entry, "encyclopedia.properties.electronic_dos") is not None
            score += functional_score.get(functional, 0)
            score += code_score.get(code_name, 0)
            if has_dos and has_bs:
                score += 10

            return (score, entry["calc_id"])

        # The calculations are first sorted by "quality"
        sorted_calc = sorted(response, key=lambda x: calc_score(x), reverse=True)

        # Get the requested representative properties
        representatives["idealized_structure"] = sorted_calc[0].calc_id
        thermo_found = False
        bs_found = False
        dos_found = False
        for calc in sorted_calc:
            if rgetattr(calc, "encyclopedia.properties.thermodynamical_properties") is not None:
                representatives["thermodynamical_properties"] = calc.calc_id
                thermo_found = True
            if rgetattr(calc, "encyclopedia.properties.electronic_band_structure") is not None:
                representatives["electronic_band_structure"] = calc.calc_id
                bs_found = True
            if rgetattr(calc, "encyclopedia.properties.electronic_dos") is not None:
                representatives["electronic_dos"] = calc.calc_id
                dos_found = True
            if thermo_found and bs_found and dos_found:
                break

        # Create result JSON
        results = []
        for entry in response:
            calc_dict = get_es_doc_values(entry, calc_prop_map)
            calc_dict["has_dos"] = calc_dict["has_dos"] is not None
            calc_dict["has_band_structure"] = calc_dict["has_band_structure"] is not None
            calc_dict["has_thermal_properties"] = calc_dict["has_thermal_properties"] is not None
            calc_dict["has_phonon_dos"] = calc_dict["has_phonon_dos"] is not None
            calc_dict["has_phonon_band_structure"] = calc_dict["has_phonon_band_structure"] is not None
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
    @api.response(200, "Metadata send", fields.Raw)
    @api.expect(statistics_query, validate=False)
    @api.marshal_with(statistics_result, skip_none=True)
    @api.doc("enc_statistics")
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
            filter=get_enc_filter() + [
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
            abort(404, message="Could not find matching calculations.")

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
        "es_source": "encyclopedia.material.idealized_structure.lattice_parameters"
    },
    "energies": {
        "es_source": "encyclopedia.properties.energies",
    },
    "mass_density": {
        "es_source": "encyclopedia.properties.mass_density",
    },
    "atomic_density": {
        "es_source": "encyclopedia.properties.atomic_density",
    },
    "cell_volume": {
        "es_source": "encyclopedia.material.idealized_structure.cell_volume"
    },
    "band_gap": {
        "es_source": "encyclopedia.properties.band_gap"
    },
    "electronic_band_structure": {
        "es_source": "encyclopedia.properties.electronic_band_structure"
    },
    "electronic_dos": {
        "es_source": "encyclopedia.properties.electronic_dos"
    },
    "phonon_band_structure": {
        "es_source": "encyclopedia.properties.phonon_band_structure"
    },
    "phonon_dos": {
        "es_source": "encyclopedia.properties.phonon_dos"
    },
    "thermodynamical_properties": {
        "es_source": "encyclopedia.properties.thermodynamical_properties"
    },
    "wyckoff_sets": {
        "arch_source": "section_metadata/encyclopedia/material/idealized_structure/wyckoff_sets"
    },
    "idealized_structure": {
        "arch_source": "section_metadata/encyclopedia/material/idealized_structure"
    },
}

calculation_property_query = api.model("calculation_query", {
    "properties": fields.List(fields.String),
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
    @api.response(200, "Metadata send", fields.Raw)
    @api.expect(calculation_property_query, validate=False)
    @api.marshal_with(calculation_property_result, skip_none=True)
    @api.doc("enc_calculation")
    def post(self, material_id, calc_id):
        """Used to return calculation details. Some properties are not
        available in the ES index and are instead read from the Archive
        directly.
        """
        # Get query parameters as json
        try:
            data = marshal(request.get_json(), calculation_property_query)
        except Exception as e:
            abort(400, message=str(e))

        s = Search(index=config.elastic.index_name)
        query = Q(
            "bool",
            filter=get_enc_filter() + [
                Q("term", encyclopedia__material__material_id=material_id),
                Q("term", calc_id=calc_id),
            ]
        )
        s = s.query(query)

        # Create dictionaries for requested properties
        references = []
        properties = data["properties"]
        arch_properties = {}
        es_properties = {}
        ref_properties = set((
            "electronic_dos",
            "electronic_band_structure",
            "thermodynamical_properties",
            "phonon_dos",
            "phonon_band_structure",
        ))
        for prop in properties:
            es_source = calculation_property_map[prop].get("es_source")
            if es_source is not None:
                es_properties[prop] = es_source
                if prop in ref_properties:
                    references.append(prop)
            arch_source = calculation_property_map[prop].get("arch_source")
            if arch_source is not None:
                arch_properties[prop] = arch_source

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
            abort(404, message="There is no material {} with calculation {}".format(material_id, calc_id))

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

        # Add results from ES
        for prop, es_source in es_properties.items():
            value = rgetattr(response[0], es_source)
            if value is not None:
                if isinstance(value, AttrDict):
                    value = value.to_dict()
                result[prop] = value

        return result, 200


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
    @api.response(204, "Report succesfully sent", fields.Raw)
    @api.expect(calculation_property_query, validate=False)
    @api.marshal_with(calculation_property_result, skip_none=True)
    @api.doc("enc_report")
    @authenticate(required=True)
    def post(self, material_id):

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
                name="webmaster", email="lauri.himanen@gmail.com", message=mail, subject='Encyclopedia error report')
        except Exception as e:
            abort(500, message="Error sending error report email.")
        print(mail)
        return "", 204


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
    upload_files = files.PublicUploadFiles(upload_id)
    with upload_files.read_archive(calc_id, access="public") as archive:
        data = archive[calc_id]
        root = EntryArchive.m_from_dict(data.to_dict())

    return root
