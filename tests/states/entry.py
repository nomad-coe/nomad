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

import json
from nomad import infrastructure, files
from nomad.utils.exampledata import ExampleData, create_entry_archive
from .archives.create_archives import archive_dft_bulk
from nomad.processing import Upload


def dft():
    """
    State containing DFT entries that can be used to e.g. test the different
    entry tabs.
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = 'dft_upload'
    entry_id = 'dft_bulk'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    entry = data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=archive_dft_bulk(),
    )

    # The archive will also be saved on disk for the GUI tests to use.
    with open('tests/states/archives/dft.json', 'w') as fout:
        json.dump(entry.m_to_dict(include_derived=True), fout, indent=2)

    data.save()


def eln():
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test').user_id
    coauthors = [infrastructure.user_management.get_user(username='scooper').user_id]
    reviewers = [infrastructure.user_management.get_user(username='ttester').user_id]
    upload = Upload(
        upload_id='eln_upload_id',
        main_author=main_author,
        coauthors=coauthors,
        reviewers=reviewers,
    )
    upload.save()
    files.StagingUploadFiles(upload_id=upload.upload_id, create=True)
    upload.staging_upload_files.add_rawfiles('examples/data/light_eln')
    upload.process_upload()
    upload.block_until_complete()


def eln_properties():
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test').user_id
    reviewers = [infrastructure.user_management.get_user(username='ttester').user_id]
    upload = Upload(
        upload_id='eln_upload_id', main_author=main_author, reviewers=reviewers
    )
    upload.save()
    files.StagingUploadFiles(upload_id=upload.upload_id, create=True)
    upload.staging_upload_files.add_rawfiles('examples/data/eln_properties')
    upload.process_upload()
    upload.block_until_complete()


def references():
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test').user_id
    coauthors = [infrastructure.user_management.get_user(username='scooper').user_id]
    reviewers = [infrastructure.user_management.get_user(username='ttester').user_id]

    upload1 = Upload(
        upload_id='references_upload_id1',
        main_author=main_author,
        coauthors=coauthors,
        reviewers=reviewers,
    )
    upload1.save()
    files.StagingUploadFiles(upload_id=upload1.upload_id, create=True)
    upload1.staging_upload_files.add_rawfiles('examples/data/references/upload1')
    upload1.process_upload()
    upload1.block_until_complete()

    upload2 = Upload(
        upload_id='references_upload_id2',
        main_author=main_author,
        coauthors=coauthors,
        reviewers=reviewers,
    )
    upload2.save()
    files.StagingUploadFiles(upload_id=upload2.upload_id, create=True)
    upload2.staging_upload_files.add_rawfiles('examples/data/references/upload2')
    upload2.process_upload()
    upload2.block_until_complete()


metadata_dict = {
    'upload_create_time': '2021-03-17T13:47:32.899000',
    'last_processing_time': '2021-03-17T15:47:32.899000',
    'nomad_version': '0.10.0',
    'nomad_commit': 'bf3c06fa',
    'comment': 'Mocked',
    'references': ['doi'],
    'main_author': {'name': 'Lauri Himanen'},
}

material_dict = {
    'material_id': 'bulk_material',
    'material_name': 'Silicon',
    'structural_type': 'bulk',
    'elements': ['Si'],
    'n_elements': 1,
    'chemical_formula_reduced': 'Si2',
    'chemical_formula_hill': 'Si2',
    'chemical_formula_anonymous': 'A2',
    'chemical_formula_descriptive': 'Si2',
    'symmetry': {
        'crystal_system': 'cubic',
        'bravais_lattice': 'cP',
        'structure_name': 'rock salt',
        'space_group_symbol': 'Fd-3m',
        'space_group_number': 227,
        'point_group': '6mm',
    },
}


def material():
    """
    Entry that contains a material.
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = 'dft_upload'
    entry_id = 'dft_bulk'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=create_entry_archive(
            metadata=metadata_dict,
            results={
                'material': material_dict,
            },
        ),
    )

    data.save()


def dos_electronic():
    """
    Entry that contains an electronic DOS.
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = 'dft_upload'
    entry_id = 'dft_bulk'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=create_entry_archive(
            metadata=metadata_dict,
            results={
                'material': material_dict,
                'properties': {
                    'available_properties': [
                        'dos_electronic',
                    ],
                    'electronic': {
                        'dos_electronic': [
                            {
                                'energies': '/run/0/calculation/0/dos_electronic/0/energies',
                                'total': [
                                    '/run/0/calculation/0/dos_electronic/0/total/0'
                                ],
                                'band_gap': [{'energy_highest_occupied': 0}],
                            }
                        ],
                    },
                },
            },
            run={
                'calculation': [
                    {
                        'dos_electronic': [
                            {
                                'energies': [0, 1e-19],
                                'total': [
                                    {
                                        'value': [0, 1e18],
                                        'normalization_factor': 1e-19,
                                        'spin': 0,
                                    }
                                ],
                                'band_gap': [
                                    {'energy_highest_occupied': 0, 'index': 0}
                                ],
                            }
                        ],
                    }
                ]
            },
        ),
    )

    data.save()


def bulk_modulus():
    """
    Entry that contains a bulk modulus.
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = 'dft_upload'
    entry_id = 'dft_bulk'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=create_entry_archive(
            metadata=metadata_dict,
            results={
                'material': material_dict,
                'properties': {
                    'available_properties': [
                        'bulk_modulus',
                    ],
                    'mechanical': {
                        'bulk_modulus': [{'type': 'murnaghan', 'value': 1}],
                    },
                },
            },
            workflow={
                'm_def': 'simulationworkflowschema.EquationOfState',
                'results': {
                    'calculation_result_ref': '/run/0/calculation/0',
                    'calculations_ref': ['/run/0/calculation/0'],
                    'equation_of_state': {
                        'energies': [0, 1],
                        'volumes': [0, 1],
                        'eos_fit': [
                            {
                                'function_name': 'murnaghan',
                                'fitted_energies': [0, 1],
                                'bulk_modulus': 1,
                            }
                        ],
                    },
                },
            },
        ),
    )

    data.save()


def trajectory():
    """
    Entry that contains a trajectory.
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = 'dft_upload'
    entry_id = 'dft_bulk'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=create_entry_archive(
            metadata=metadata_dict,
            results={
                'material': material_dict,
                'properties': {
                    'available_properties': ['trajectory'],
                    'thermodynamic': {
                        'trajectory': [
                            {
                                'pressure': {
                                    'value': [0],
                                    'time': [0],
                                },
                                'temperature': {
                                    'value': [0],
                                    'time': [0],
                                },
                                'volume': {
                                    'value': [0],
                                    'time': [0],
                                },
                                'available_properties': [
                                    'pressure',
                                    'temperature',
                                    'volume',
                                ],
                                'methodology': {
                                    'molecular_dynamics': {
                                        'time_step': 1e-15,
                                        'ensemble_type': 'NVT',
                                    }
                                },
                            }
                        ]
                    },
                },
            },
            run={
                'calculation': [
                    {
                        'pressure': 0,
                        'temperature': 0,
                        'volume': 0,
                    }
                ]
            },
            workflow={
                'm_def': 'simulationworkflowschema.MolecularDynamics',
                'results': {
                    'calculation_result_ref': '/run/0/calculation/0',
                    'calculations_ref': ['/run/0/calculation/0'],
                },
                'method': {
                    'integration_timestep': 1e-15,
                    'thermodynamic_ensemble': 'NVT',
                },
            },
        ),
    )

    data.save()


def dos_phonon():
    """
    Entry that contains a phonon DOS.
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = 'dft_upload'
    entry_id = 'dft_bulk'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=create_entry_archive(
            metadata=metadata_dict,
            results={
                'material': material_dict,
                'properties': {
                    'available_properties': [
                        'dos_phonon',
                    ],
                    'vibrational': {
                        'dos_phonon': {
                            'energies': '/run/0/calculation/0/dos_phonon/0/energies',
                            'total': ['/run/0/calculation/0/dos_phonon/0/total/0'],
                        },
                    },
                },
            },
            run={
                'calculation': [
                    {
                        'dos_phonon': [
                            {
                                'energies': [0, 1e-19],
                                'total': [
                                    {'value': [0, 1e18], 'normalization_factor': 1e-19}
                                ],
                            }
                        ],
                    }
                ]
            },
            workflow={
                'm_def': 'simulationworkflowschema.Phonon',
                'results': {
                    'calculation_result_ref': '/run/0/calculation/0',
                    'calculations_ref': ['/run/0/calculation/0'],
                    'thermodynamics': {
                        'heat_capacity_c_v': [0, 1],
                        'vibrational_free_energy_at_constant_volume': [0, 1],
                        'temperature': [0, 100],
                    },
                },
            },
        ),
    )

    data.save()


def rdf():
    """
    Entry that contains an RDF.
    """
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = 'dft_upload'
    entry_id = 'dft_bulk'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=create_entry_archive(
            metadata=metadata_dict,
            results={
                'material': material_dict,
                'properties': {
                    'available_properties': [
                        'radial_distribution_function',
                    ],
                    'structural': {
                        'radial_distribution_function': [
                            {
                                'type': 'molecular',
                                'label': '0-0',
                                'bins': [
                                    8.408812522888184e-12,
                                ],
                                'n_bins': 1,
                                'value': [
                                    0.0,
                                ],
                                'frame_start': 0,
                                'frame_end': 40,
                                'methodology': {
                                    'molecular_dynamics': {
                                        'time_step': 2.5e-16,
                                        'ensemble_type': 'NVT',
                                    }
                                },
                            }
                        ]
                    },
                },
            },
        ),
    )

    data.save()


def plotly():
    infrastructure.setup()
    main_author = infrastructure.user_management.get_user(username='test').user_id
    upload = Upload(upload_id='plotly_upload_id', main_author=main_author)
    upload.save()
    files.StagingUploadFiles(upload_id=upload.upload_id, create=True)
    upload.staging_upload_files.add_rawfiles('examples/data/plotly')
    upload.process_upload()
    upload.block_until_complete()
