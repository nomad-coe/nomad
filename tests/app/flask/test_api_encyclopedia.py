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
import time
import pytest
import click.testing
import json
import elasticsearch

from nomad import config
from nomad.cli import cli
from nomad import processing as proc, infrastructure

from tests.app.flask.test_app import BlueprintClient

silicon_id = "fh3UBjhUVm4nxzeRd2JJuqw5oXYa"


def exists(value, expected_type=str, is_not=set(["", None])):
    return type(value) == expected_type and value not in is_not


def validate_material(material):
    assert exists(material["material_id"])
    assert exists(material["formula"])
    assert exists(material["formula_reduced"])
    material_type = material["material_type"]
    assert exists(material_type)
    if material_type == "bulk":
        assert exists(material["material_name"])
        assert exists(material["has_free_wyckoff_parameters"], bool)
        assert exists(material["strukturbericht_designation"])
        assert exists(material["bravais_lattice"])
        assert exists(material["crystal_system"])
        assert exists(material["point_group"])
        assert exists(material["space_group_number"], int)
        assert exists(material["space_group_international_short_symbol"])
        assert exists(material["structure_type"])
        assert exists(material["structure_prototype"])


@pytest.fixture(scope='function')
def api(client):
    return BlueprintClient(client, '/api/encyclopedia')


def upload(filepath, publish, test_user_bravado_client, proc_infra, test_user_auth, api):

    # Perform a test upload
    arguments = ['client', 'upload', '--offline', '--name', filepath, filepath]
    if publish:
        arguments.append("--publish")
    click.testing.CliRunner().invoke(
        cli,
        arguments,
        catch_exceptions=False
    )
    upload = proc.Upload.objects(name=filepath).first()
    upload_id = upload["upload_id"]
    upload.block_until_complete(interval=0.2)

    return upload_id


@pytest.fixture(scope='function')
def enc_upload(test_user_bravado_client, proc_infra, test_user_auth, api, mongo_infra):
    upload('tests/data/api/enc_public.zip', True, test_user_bravado_client, proc_infra, test_user_auth, api)
    upload('tests/data/api/enc_private_material.zip', False, test_user_bravado_client, proc_infra, test_user_auth, api)
    upload('tests/data/api/enc_si_private.zip', False, test_user_bravado_client, proc_infra, test_user_auth, api)
    upload_id = upload('tests/data/api/enc_si_embargo.zip', True, test_user_bravado_client, proc_infra, test_user_auth, api)

    # Place upload entries on embargo in MongoDB
    calculations = mongo_infra["test_db"]['calc']
    calculations.update_many(
        {'upload_id': upload_id},
        {'$set': {'metadata.with_embargo': True}}
    )

    # Place upload entries on embargo in ES
    embargoed = proc.Upload.get(upload_id)
    with embargoed.entries_metadata() as calcs:
        def elastic_updates():
            for calc in calcs:
                entry = calc.a_elastic.create_index_entry()
                entry.with_embargo = True
                entry = entry.to_dict(include_meta=True)
                source = entry.pop('_source')
                entry['doc'] = source
                entry['_op_type'] = 'update'
                entry['_index'] = 'nomad_fairdi_calcs_test'
                yield entry

        elasticsearch.helpers.bulk(infrastructure.elastic_client, elastic_updates())
        # The indices need to be refreshed after the update
        infrastructure.elastic_client.indices.refresh(config.elastic.index_name)

    # Index materials
    click.testing.CliRunner().invoke(
        cli,
        ['admin', 'index-materials', "--source=es"],
        catch_exceptions=False)

    # Small wait time in order for ES indexing to finish up
    time.sleep(1)


class TestEncyclopedia():

    # @pytest.mark.skip(reason='this still fails due to metainfo refactor and needs fixing')
    def test_material(self, enc_upload, elastic_infra, api, test_user_auth):
        # Correctly found material returns all required values.
        rv = api.get('/materials/{}'.format(silicon_id))
        assert rv.status_code == 200
        material = rv.json
        validate_material(material)

        # Missing material causes 404.
        rv = api.get('/materials/does_not_exist')
        assert rv.status_code == 404

        # Empty search should return all visible materials
        rv = api.post(
            '/materials/',
            data="{}",
            content_type='application/json'
        )
        assert rv.status_code == 200
        results = rv.json['results']
        assert len(results) == 2
        for material in results:
            validate_material(material)

        # Test formula search
        rv = api.post(
            '/materials/',
            data=json.dumps({"search_by": {"formula": "Si"}}),
            content_type='application/json'
        )
        assert rv.status_code == 200
        rv = api.post(
            '/materials/',
            data=json.dumps({"search_by": {"formula": "mgbi2"}}),
            content_type='application/json'
        )
        assert rv.status_code == 400
        rv = api.post(
            '/materials/',
            data=json.dumps({"search_by": {"formula": "mgBi2"}}),
            content_type='application/json'
        )
        assert rv.status_code == 400
        rv = api.post(
            '/materials/',
            data=json.dumps({"search_by": {"formula": "Mgbi2"}}),
            content_type='application/json'
        )
        assert rv.status_code == 400
        rv = api.post(
            '/materials/',
            data=json.dumps({"search_by": {"formula": "MgaBi2"}}),
            content_type='application/json'
        )
        assert rv.status_code == 400

        # Test that searches across calculations work as intended
        rv = api.post(
            '/materials/',
            data=json.dumps({
                "has_band_structure": True,
                "has_dos": True,
            }),
            content_type='application/json'
        )
        results = rv.json['results']
        assert len(results) == 1
        rv = api.post(
            '/materials/',
            data=json.dumps({
                "code_name": ["Quantum Espresso"],
                "functional_type": ["GGA"],
                "has_band_structure": True,
            }),
            content_type='application/json'
        )
        results = rv.json['results']
        assert len(results) == 1

        # Test element search
        rv = api.post(
            '/materials/',
            data=json.dumps({"search_by": {"elements": ["Si"]}}),
            content_type='application/json'
        )
        assert rv.status_code == 200

        # Test that searches within calculations work as intended
        rv = api.post(
            '/materials/',
            data=json.dumps({
                "search_by": {
                    "restricted": True,
                },
                "has_band_structure": True,
                "has_dos": True,
            }),
            content_type='application/json'
        )
        results = rv.json['results']
        assert len(results) == 0
        rv = api.post(
            '/materials/',
            data=json.dumps({
                "search_by": {
                    "restricted": True,
                },
                "code_name": ["Quantum Espresso"],
                "functional_type": ["GGA"],
                "has_band_structure": True,
            }),
            content_type='application/json'
        )
        results = rv.json['results']
        assert len(results) == 0

        # Test that EOS groups are found.
        rv = api.get('/materials/{}/groups'.format(silicon_id))
        assert rv.status_code == 200
        groups = rv.json
        groups_eos = groups['groups_eos']
        assert len(groups_eos) == 1
        eos_id, group = list(groups_eos.items())[0]
        exists(eos_id)
        assert len(group) == 5
        for calc_id in group:
            exists(calc_id)

        # Test query for a specific EOS group.
        rv = api.get('/materials/{}/groups/eos/{}'.format(silicon_id, eos_id))
        assert rv.status_code == 200
        group = rv.json
        assert len(group['calculations']) == 5
        assert len(group['energies']) == 5
        assert len(group['volumes']) == 5

        # TODO: parameter variation groups should be tested at some point.
        # Test that parameter variation groups are found.
        # rv = api.get('/materials/{}/groups'.format(silicon_id))
        # assert rv.status_code == 200
        # groups = rv.json
        # groups_par = groups["groups_par"]
        # assert len(groups_par) == 1
        # par_id, group = list(groups_par.items())[0]
        # exists(par_id)
        # assert len(group) == 2
        # for calc_id in group:
        # exists(calc_id)

        # rv = api.get('/materials/{}/groups/par/{}'.format(silicon_id, par_id))
        # assert rv.status_code == 200
        # group = rv.json
        # assert len(group['calculations']) == 2
        # assert len(group['energies']) == 2
        # assert len(group['volumes']) == 2

        # Test suggestions
        rv = api.get('/suggestions?property=structure_type')
        assert rv.status_code == 200
        structure_types = rv.json
        assert structure_types["structure_type"] == ["diamond"]
        rv = api.get('/suggestions?property=code_name')
        assert rv.status_code == 200
        code_names = rv.json
        assert set(code_names["code_name"]) == set(["exciting", "Quantum Espresso", "VASP", "FHI-aims"])

        # Test calculations: embargoed and private entries should not show up here
        rv = api.get('/materials/{}/calculations'.format(silicon_id))
        assert rv.status_code == 200
        calculations = rv.json
        assert calculations["total_results"] == 9
        assert len(calculations["results"]) == calculations["total_results"]
        calc_ids = [x["calc_id"] for x in calculations["results"]]
        for calc in calculations["results"]:
            assert exists(calc["calc_id"])
            assert exists(calc["upload_id"])
            assert exists(calc["code_name"])
            assert exists(calc["code_version"])
            assert exists(calc["functional_type"])
            assert exists(calc["basis_set_type"])
            assert exists(calc["core_electron_treatment"])
            assert exists(calc["run_type"])
            assert exists(calc["has_dos"], bool)
            assert exists(calc["has_band_structure"], bool)
            assert exists(calc["has_thermal_properties"], bool)

        # Test statistics
        rv = api.post(
            '/materials/{}/statistics'.format(silicon_id),
            data=json.dumps({
                "calculations": calc_ids,
                "properties": [
                    "cell_volume",
                    "atomic_density",
                    "mass_density",
                    "lattice_a",
                    "lattice_b",
                    "lattice_c",
                    "alpha",
                    "beta",
                    "gamma",
                    "band_gap",
                ],
                "n_histogram_bins": 3,
            }),
            content_type='application/json'
        )
        assert rv.status_code == 200

        # Test fetching information about a specific calculation
        rv = api.post(
            '/materials/{}/calculations/{}'.format(silicon_id, calc_ids[0]),
            data=json.dumps({"properties": [
                "lattice_parameters",
                "energies",
                "mass_density",
                "atomic_density",
                "cell_volume",
                "wyckoff_sets",
                "idealized_structure",
                "band_gap",
                "electronic_band_structure",
                "electronic_dos",
                "phonon_band_structure",
                "phonon_dos",
                "thermodynamical_properties",
            ]}),
            content_type='application/json'
        )
        assert rv.status_code == 200
        calc = rv.json

        # Test that completely private materials only become visible after
        # authentication
        rv = api.post(
            '/materials/',
            data=json.dumps({"search_by": {"elements": ["B"]}}),
            content_type='application/json',
        )
        results = rv.json['results']
        assert len(results) == 0
        rv = api.post(
            '/materials/',
            data=json.dumps({"search_by": {"elements": ["B"]}}),
            content_type='application/json',
            headers=test_user_auth,
        )
        results = rv.json['results']
        assert len(results) == 1
        private_material = results[0]
        private_material_id = private_material["material_id"]
        rv = api.get('/materials/{}/calculations'.format(private_material_id), headers=test_user_auth)
        private_calc_ids = [x["calc_id"] for x in rv.json["results"]]
        for headers, code in [({}, 404), (test_user_auth, 200)]:
            rv = api.get('/materials/{}'.format(private_material_id), headers=headers)
            assert rv.status_code == code
            rv = api.get('/materials/{}/groups'.format(private_material_id), headers=headers)
            assert rv.status_code == code
            rv = api.get('/materials/{}/calculations'.format(private_material_id), headers=headers)
            assert rv.status_code == code
            rv = api.post(
                '/materials/{}/statistics'.format(private_material_id),
                data=json.dumps({
                    "calculations": private_calc_ids,
                    "properties": [
                        "cell_volume",
                        "atomic_density",
                        "mass_density",
                        "lattice_a",
                        "lattice_b",
                        "lattice_c",
                        "alpha",
                        "beta",
                        "gamma",
                        "band_gap",
                    ],
                    "n_histogram_bins": 3,
                }),
                content_type='application/json',
                headers=headers,
            )
            assert rv.status_code == code

        # Test that invalid query parameters raise code 400

    def test_complex_search(self, enc_upload, elastic_infra, api, test_user_auth):
        # Test an elaborate boolean query for elements
        query = json.dumps({"query": """(
          ( elements HAS ALL "Si", "O" OR elements HAS ALL "Ge", "O" ) OR
          ( elements HAS ALL "Si", "N" OR elements HAS ALL "Ge", "N" )
        )"""})
        rv = api.post('/materials/', data=query, content_type='application/json')
        assert rv.status_code == 200
        results = rv.json['results']
        assert len(results) == 0

        # Test that there are no issues with the custom Optimade grammar
        # containing boolean values. See discussion at
        # https://github.com/Materials-Consortia/OPTIMADE/issues/345
        query = json.dumps({"query": 'has_band_structure=TRUE'})
        rv = api.post('/materials/', data=query, content_type='application/json')
        assert rv.status_code == 200
        results = rv.json['results']
        assert len(results) == 1

        # Test that completely private materials only become visible after
        # authentication
        query = json.dumps({"query": 'elements HAS ALL "B"'})
        rv = api.post('/materials/', data=query, content_type='application/json')
        results = rv.json['results']
        assert rv.status_code == 200
        assert len(results) == 0
        rv = api.post('/materials/', data=query, content_type='application/json', headers=test_user_auth)
        assert rv.status_code == 200
        results = rv.json['results']
        assert len(results) == 1
