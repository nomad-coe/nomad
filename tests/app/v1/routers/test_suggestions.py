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
These are the tests for all API operations below ``suggestions``. The tests are organized
using the following type of methods: fixtures, ``perform_*_test``, ``assert_*``, and
``test_*``. While some ``test_*`` methods test individual API operations, some
test methods will test multiple API operations that use common aspects like
supporting queries, pagination, or the owner parameter. The test methods will use
``perform_*_test`` methods as an parameter. Similarly, the ``assert_*`` methods allow
to assert for certain aspects in the responses.
"""

import pytest
from tests.utils import ExampleData
from nomad.metainfo.elasticsearch_extension import entry_type
from .common import assert_response


@pytest.fixture(scope="module")
def example_data_suggestions(elastic_module, raw_files_module, mongo_module, test_user, other_test_user, normalized):
    data = ExampleData(main_author=test_user)
    upload_id = "suggestions_upload"

    data.create_upload(
        upload_id=upload_id,
        published=True
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id="suggestions_entry_1",
        mainfile="test_content/test_entry/mainfile.json",
        results={
            "material": {
                "material_name": "alpha/beta",
                "elements": ["C", "H", "Br"],
                "chemical_formula_hill": "C2H5Br",
                "chemical_formula_anonymous": "A2B5C",
                "chemical_formula_descriptive": "C2H5Br",
                "chemical_formula_reduced": "C2H5Br",
                # A duplicate value is given to test that only unique values
                # are suggested.
                "functional_type": ["semimetal", "semiconductor", "semimetal"],
                "symmetry": {
                    "crystal_system": "cubic",
                    "structure_name": "rock salt",
                    "prototype_aflow_id": "AB_cF8_225_a_b",
                }
            },
            "method": {
                "simulation": {
                    "program_name": "test_name",
                    "program_version": "10.12",
                    "dft": {
                        "xc_functional_names": ["GGA_X_PBE_SOL", "GGA_C_PBE_SOL"]
                    }
                }
            },
            "properties": {
                "mechanical": {
                    "bulk_modulus": [
                        {
                            "type": "birch_euler",
                            "value": 1,
                        }
                    ]
                }
            }
        }
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id="suggestions_entry_2",
        mainfile="test_content/test_entry/mainfile.json",
        results={
            "material": {
                "chemical_formula_hill": "ClNa",
            },
        }
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id="suggestions_entry_3",
        mainfile="test_content/test_entry/mainfile.json",
        results={
            "material": {
                "chemical_formula_hill": "Ni2O2",
            },
        }
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id="suggestions_entry_4",
        mainfile="test_content/test_entry/mainfile.json",
        results={
            "material": {
                "chemical_formula_hill": "Mg2O2",
            },
        }
    )

    data.save()

    yield

    data.delete()
    from nomad.search import search
    assert search(query=dict(upload_id=upload_id)).pagination.total == 0


def run_query(quantities, input, client):
    body = {
        "input": input,
        "quantities": [quantities] if isinstance(quantities, str) else quantities
    }
    response = client.post("suggestions", json=body, headers={})
    return response


def assert_suggestions(quantity, input, output, client):
    if isinstance(input, str):
        input = [input]
    for x in input:
        response = run_query(quantity, x, client)
        assert_response(response, 200)

        response = response.json()
        suggestions = response[quantity]
        suggestion_values = set([suggestion["value"] for suggestion in suggestions])
        if isinstance(output, str):
            output = [output]
        output = set(output)
        assert output.issubset(suggestion_values)
        assert len(suggestions) == len(output)


def test_suggestions_unknown(client, example_data_suggestions):
    """Tests that trying to get suggestions for an unregistered quantity gives
    the correct status code.
    """
    response = run_query("does.not.exist", "cu", client)
    assert_response(response, 422)


def test_suggestions_all(client, example_data_suggestions):
    """Test that running the query against all defined suggestion values works
    correctly.
    """
    response = run_query(list(entry_type.suggestions), "cu", client)
    assert_response(response, 200)


@pytest.mark.parametrize("quantity, input, output", [
    # "simple" tokenizer
    ("results.material.symmetry.crystal_system", "cu", "cubic"),

    # "default" tokenizer
    ("results.material.symmetry.structure_name", ["sa", "ro"], "rock salt"),  # Whitespace tokenization
    ("results.method.simulation.program_name", ["te", "na"], "test_name"),  # Underscore tokenization
    ("results.material.material_name", ["al", "be"], "alpha/beta"),  # Slash tokenization
    ("results.method.simulation.program_version", ["10", "12"], "10.12"),  # Dot tokenization and numbers
    # Input that spans across several tokenized words and does not start from the beginning
    ("results.method.simulation.dft.xc_functional_names", "PBE_SOL", ["GGA_C_PBE_SOL", "GGA_X_PBE_SOL"]),

    # Only "formula" tokenizer
    ("results.material.chemical_formula_anonymous", ["A2", "B5", "C"], "A2B5C"),
    ("results.material.chemical_formula_descriptive", ["C2", "H5", "Br"], "C2H5Br"),
    ("results.material.chemical_formula_reduced", ["C2", "H5", "Br"], "C2H5Br"),

    # "formula" tokenizer and "formula" variants. Note that the best matching
    # variant is returned
    ("results.material.chemical_formula_hill", ["Na", "NaC", "NaCl"], "NaCl"),
    ("results.material.chemical_formula_hill", ["Cl", "ClN", "ClNa"], "ClNa"),

    # Tests that all matches are returned even if they share the same matched token
    ("results.material.chemical_formula_hill", "O2", ["Ni2O2", "Mg2O2"]),

    # Nested fields
    ("results.properties.mechanical.bulk_modulus.type", "euler", "birch_euler"),

    # Fields with multiple values
    ("results.material.elements", "Br", "Br"),
    ("results.material.functional_type", "semic", "semiconductor"),
    ("results.material.functional_type", "semim", "semimetal"),
    ("results.material.functional_type", "semi", ["semiconductor", "semimetal"]),
])
def test_suggestions_quantities(quantity, input, output, client, example_data_suggestions):
    assert_suggestions(quantity, input, output, client)
