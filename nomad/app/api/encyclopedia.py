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

'''
The encyclopedia API of the nomad@FAIRDI APIs.
'''
import re
from math import gcd
from functools import reduce

import numpy as np
from flask_restplus import Resource, abort, fields, marshal
from flask import request
from elasticsearch_dsl import Search, Q

from nomad import config
from nomad.units import ureg
from nomad.atomutils import get_hill_decomposition
from .api import api

ns = api.namespace('encyclopedia', description='Access encyclopedia metadata.')
re_formula = re.compile(r"([A-Z][a-z]?)(\d*)")


def add_result(result, key, function, default=""):
    """Convenience function that attempts to add a value from the ElasticSearch
    result into the given result object. Upon failing returns the specified
    default value.
    """
    try:
        value = function()
    except Exception:
        value = default
    result[key] = value


def get_material(es_doc):
    """Used to form a material definition from the given ElasticSearch root
    document.
    """
    result = {}
    add_result(result, "material_id", lambda: es_doc.encyclopedia.material.material_id, ""),
    add_result(result, "bravais_lattice", lambda: es_doc.encyclopedia.material.bulk.bravais_lattice, ""),
    add_result(result, "crystal_system", lambda: es_doc.encyclopedia.material.bulk.crystal_system, "")
    add_result(result, "formula", lambda: es_doc.encyclopedia.material.formula, "")
    add_result(result, "formula_reduced", lambda: es_doc.encyclopedia.material.formula_reduced, "")
    add_result(result, "material_name", lambda: es_doc.encyclopedia.material.material_name, "")
    add_result(result, "point_group", lambda: es_doc.encyclopedia.material.bulk.point_group, "")
    add_result(result, "space_group", lambda: es_doc.encyclopedia.material.bulk.space_group_number, "")
    add_result(result, "structure_type", lambda: es_doc.encyclopedia.material.bulk.structure_type, "")
    add_result(result, "system_type", lambda: es_doc.encyclopedia.material.material_type, "")

    return result


material_query = api.parser()
material_query.add_argument('material_id', type=str, help='Identifier for the searched material.', location='args')
material_result = api.model('material_result', {
    "bravais_lattice": fields.String,
    "crystal_system": fields.String,
    "formula": fields.String,
    "formula_reduced": fields.String,
    "material_name": fields.String,
    "point_group": fields.String,
    "space_group": fields.Integer(),
    "structure_type": fields.String,
    "system_type": fields.String,
})


@ns.route('/materials/<string:material_id>')
class EncMaterialResource(Resource):
    @api.response(404, 'The material does not exist')
    @api.response(200, 'Metadata send', fields.Raw)
    @api.doc('material/<material_id>')
    @api.expect(material_query)
    @api.marshal_with(material_result)
    def get(self, material_id):
        """Used to retrive basic information related to the specified material.
        """

        # Find the first public entry with this material id and take
        # information from there. In principle all other entries should have
        # the same information.
        s = Search(index=config.elastic.index_name)

        # Since we are looking for an exact match, we use filter context
        # together with term search for speed (instead of query context and
        # match search)
        query = Q(
            'bool',
            filter=[
                Q('term', published=True),
                Q('term', with_embargo=False),
                Q('term', encyclopedia__material__material_id=material_id),
            ]
        )
        s = s.query(query)
        response = s.execute()

        # No such material
        if len(response) == 0:
            abort(404, message='There is no material {}'.format(material_id))

        # Create result JSON
        entry = response[0]
        result = get_material(entry)

        return result, 200


range_query = api.model('range_query', {
    "max": fields.Float,
    "min": fields.Float,
})
materials_query = api.model('materials_input', {
    'search_by': fields.Nested(api.model('search_query', {
        "exclusive": fields.Boolean(default=False),
        "formula": fields.String,
        "element": fields.String,
        "page": fields.Integer(default=1),
        "per_page": fields.Integer(default=25),
        "pagination": fields.Boolean,
    })),
    'material_name': fields.List(fields.String),
    'structure_type': fields.List(fields.String),
    'space_group': fields.List(fields.Integer),
    'system_type': fields.List(fields.String),
    'crystal_system': fields.List(fields.String),
    'band_gap': fields.Nested(range_query, description="Band gap range in eV."),
    'band_gap_direct': fields.Boolean,
    'has_band_structure': fields.Boolean,
    'has_dos': fields.Boolean,
    'has_fermi_surface': fields.Boolean,
    'has_thermal_properties': fields.Boolean,
    'functional_type': fields.List(fields.String),
    'basis_set_type': fields.List(fields.String),
    'code_name': fields.List(fields.String),
    'mass_density': fields.Nested(range_query, description="Mass density range in kg / m ** 3."),
})
materials_result = api.model('materials_result', {
    'pages': fields.Integer(required=True),
    'results': fields.List(fields.Nested(material_result)),
    'total_results': fields.Integer(allow_null=False),
})


@ns.route('/materials')
class EncMaterialsResource(Resource):
    @api.response(404, 'No materials found')
    @api.response(400, 'Bad request')
    @api.response(200, 'Metadata send', fields.Raw)
    @api.expect(materials_query, validate=False)
    @api.marshal_with(materials_result)
    @api.doc('materials')
    def post(self):
        """Used to query a list of materials with the given search options.
        """
        # Get query parameters as json
        try:
            data = marshal(request.get_json(), materials_query)
        except Exception as e:
            abort(400, message=str(e))

        s = Search(index=config.elastic.index_name)
        filters = []
        must_nots = []
        musts = []

        # Add term filters
        filters.append(Q('term', published=True))
        filters.append(Q('term', with_embargo=False))

        def add_terms_filter(source, target, query_type="terms"):
            if data[source]:
                filters.append(Q(query_type, **{target: data[source]}))

        add_terms_filter("material_name", "encyclopedia.material.material_name")
        add_terms_filter("structure_type", "encyclopedia.material.bulk.structure_type")
        add_terms_filter("space_group", "encyclopedia.material.bulk.space_group_number")
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
            names, counts = get_hill_decomposition(element_list)
            greatest_common_divisor = reduce(gcd, counts)
            reduced_counts = np.array(counts) / greatest_common_divisor
            query_string = []
            for name, count in zip(names, reduced_counts):
                if count == 1:
                    query_string.append(name)
                else:
                    query_string.append("{}{}".format(name, count))
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

        # Prepare the final boolean query that combines the different queries
        filter_query = Q('bool', filter=filters, must_not=must_nots, must=musts)
        s = s.query(filter_query)

        # Execute query
        response = s.execute()

        # No matches
        if len(response) == 0:
            abort(404, message='No materials found for the given search criteria.')

        # Create final result dictionary
        result_list = [get_material(es_doc) for es_doc in response]
        result = {
            "total_results": len(result_list),
            "pages": None,
            "results": result_list,
        }
        return result, 200

# @ns.route('/esmaterials')
# class EncESMaterialsResource(Resource):
    # @api.response(404, 'No materials found')
    # @api.response(200, 'Metadata send', fields.Raw)
    # @api.doc('materials')
    # def post(self):
        # """Used to query a list of materials with the given ElasticSearch JSON
        # query.
        # """
