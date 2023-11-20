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
from datetime import datetime, timedelta
from nomad import infrastructure
from nomad.utils import create_uuid
from nomad.units import ureg
from nomad.atomutils import chemical_symbols
from nomad.utils.exampledata import ExampleData

material_h2o = {
    'structural_type': 'molecule / cluster',
    'elements': ['C', 'H'],
    'chemical_formula_hill': 'CH3',
    'chemical_formula_anonymous': 'AB2',
    'chemical_formula_descriptive': 'CH3',
    'chemical_formula_reduced': 'CH3',
}
material_high_entropy_alloy = {
    'material_id': 'aLzu5XmhATsTfEA4sGIiL45caCMq',
    'structural_type': 'bulk',
    'elements': ['Hf', 'Nb', 'Ta', 'Ti', 'Zr'],
    'chemical_formula_hill': 'HfNbTaTiZr',
    'chemical_formula_anonymous': 'ABCDE',
    'chemical_formula_descriptive': 'HfNbTaTiZr',
    'chemical_formula_reduced': 'HfNbTaTiZr',
}
material_graphene = {
    'material_id': 'jPcWsq-Rb0gtgx2krJDmLvQUxpcL',
    'structural_type': '2D',
    'elements': ['C'],
    'chemical_formula_hill': 'C',
    'chemical_formula_anonymous': 'A',
    'chemical_formula_descriptive': 'C',
    'chemical_formula_reduced': 'C',
}
material_perovskite = {
    "material_name": "perovskite",
    "structural_type": "bulk",
    "functional_type": ["semiconductor", "solar cell"],
    "elements": ["C", "H", "I", "N", "Pb"],
    "chemical_formula_descriptive": "MAPbI3",
    "chemical_formula_reduced": "H6Pb1C1I3N1",
    "chemical_formula_hill": "CH6I3NPb"
}
material_movo = {
    'material_id': 'mock-material-id',
    "material_name": "movo",
    "elements": ["V", "Mo", "O"],
    "chemical_formula_descriptive": "MoVO",
    "chemical_formula_reduced": "MoVO",
    "chemical_formula_hill": "MoVO"
}
material_mof = {
    'material_id': 'mock-material-id',
    "material_name": "movo",
    "elements": ["Fe", "C", "O"],
    "chemical_formula_descriptive": "COFe",
    "chemical_formula_reduced": "COFe",
    "chemical_formula_hill": "COFe",
    "topology": [{'label': 'MOF'}]
}
method_dft_vasp = {
    'simulation': {
        'program_name': 'VASP',
        'dft': {
            'xc_functional_names': ['GGA_X_PBE_SOL', 'GGA_C_PBE_SOL']
        }
    }
}
method_dft_exciting = {
    'simulation': {
        'program_name': 'exciting',
        'dft': {
            'xc_functional_names': ['LDA_X_PZ', 'LDA_C_PZ']
        }
    }
}
method_eels = {
    'method_name': 'EELS'
}
properties_solar_cell = {
    "available_properties": ["solar_cell", "electronic.band_structure_electronic.band_gap"],
    'electronic': {
        'band_structure_electronic': [{
            'band_gap': [{
                'value': 2.5634826144e-19
            }]
        }]
    },
    "optoelectronic": {
        "solar_cell": {
            "efficiency": 15.9,
            "fill_factor": 0.78,
            "open_circuit_voltage": 1.0,
            "short_circuit_current_density": 203.0,
            "illumination_intensity": 1000.0,
            "device_area": 1.25e-05,
            "device_architecture": "pin",
            "device_stack": ["SLG", "ITO", "PEDOT:PSS", "Perovskite", "PCBM-60", "BCP", "Ag"],
            "absorber": ["MAPbI"],
            "absorber_fabrication": ["Spin-coating"],
            "electron_transport_layer": ["PCBM-60", "BCP"],
            "hole_transport_layer": ["PEDOT:PSS"],
            "substrate": ["SLG", "ITO"],
            "back_contact": ["Ag"]
        }
    }
}
properties_catalysis = {
    "catalytic": {
        "catalyst_characterization": {
            "surface_area": 65.56,
            "method_surface_area": "BET"
        },
        "catalyst_synthesis": {
            "catalyst_type": "bulk catalyst",
            "preparation_method": "hydrothermal"
        }
    }
}
eln_catalysis = {
    "sections": [
        "CatalystSample"
    ],
    "names": [
        "MoVOx sample"
    ]
}


def search():
    '''
    Used to construct an API state that is suitable for several different kinds of
    search tests. Constructs an upload with multiple entries containing a wide
    variety of data.
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = create_uuid()
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)

    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_h2o,
            'method': method_dft_vasp,
            'properties': {}
        },
        quantities=['results.method.simulation.program_name']
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_graphene,
            'method': method_dft_exciting,
            'properties': {}
        },
        quantities=['results.method.simulation.program_name']
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_high_entropy_alloy,
            'method': method_dft_vasp,
            'properties': {}
        },
        quantities=['results.method.simulation.program_name']
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_perovskite,
            'method': {},
            'properties': properties_solar_cell
        },
        sections=['nomad.datamodel.results.SolarCell'],
        quantities=['data']
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/movo.json',
        results={
            'material': material_movo,
            'properties': properties_catalysis,
            'eln': eln_catalysis
        },
        quantities=['results.properties.catalytic']
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_perovskite,
            'method': method_eels,
            'properties': {}
        },
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_mof,
            'method': method_dft_vasp,
            'properties': {}
        },
    )

    data.save()


def histograms():
    '''
    Used to construct an API state that is suitable for testing histograms.
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = create_uuid()
    upload_create_time = datetime.utcfromtimestamp(1585872000)  # 3/4/2020 00:00 GMT
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0, upload_create_time=upload_create_time)

    # Create entries with different number of elements for testing the
    # discretized n_elements field
    for i in [1, 2, 2, 3, 4, 5]:
        data.create_entry(
            upload_id=upload_id,
            entry_id=create_uuid(),
            mainfile='upload/archive.json',
            results={
                'material': {'elements': chemical_symbols(range(1, i + 1))},
                'method': {},
                'properties': {}
            }
        )

    # Create entries with different band gaps for testing continuous histograms
    for i in [0, 0.1, 0.5, 1, 2.0]:
        data.create_entry(
            upload_id=upload_id,
            entry_id=create_uuid(),
            mainfile='upload/archive.json',
            results={
                'material': material_h2o,
                'method': {},
                'properties': {
                    'electronic': {
                        'band_structure_electronic': [{
                            'band_gap': [{
                                'value': (i * ureg.eV).to(ureg.joule).magnitude
                            }]
                        }]
                    }
                }
            }
        )

    # Create dummy uploads to test the timestamp histogram
    for i in range(1, 5):
        upload_id = create_uuid()
        upload_create_time += timedelta(minutes=1)
        data.create_upload(upload_id=upload_id, published=True, embargo_length=0, upload_create_time=upload_create_time)
        data.create_entry(
            upload_id=upload_id,
            entry_id=create_uuid(),
            mainfile='upload/archive.json',
            results={
                'material': {'elements': ['H']},
                'method': {},
                'properties': {}
            }
        )

    data.save()


def histograms_one_value():
    '''
    Creates a state with one value for each type of histogram.
    '''
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = create_uuid()
    upload_create_time = datetime.utcfromtimestamp(1585872000)  # 3/4/2020 00:00 GMT
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0, upload_create_time=upload_create_time)

    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': {'elements': ['H']},
            'method': {},
            'properties': {
                'electronic': {
                    'band_structure_electronic': [{
                        'band_gap': [{
                            'value': (0.5 * ureg.eV).to(ureg.joule).magnitude
                        }]
                    }]
                }
            }
        }
    )

    data.save()
