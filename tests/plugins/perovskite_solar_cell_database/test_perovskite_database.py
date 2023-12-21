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
import os

from nomad.datamodel.datamodel import EntryArchive
from nomad.datamodel.context import ClientContext
from nomad.utils.exampledata import ExampleData


def test_perovskite_solar_cell_plugin_processing(
    raw_files_function, no_warn, test_user, mongo_function
):
    directory = 'tests/data/plugins/perovskite_solar_cell_database'
    mainfile = 'example.archive.json'
    upload_id = 'test_upload_id'
    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id=upload_id, published=False)
    context = ClientContext(local_dir=directory, upload_id=upload_id)

    test_archive = data.create_entry_from_file(
        upload_id=upload_id,
        mainfile=os.path.join(directory, mainfile),
        entry_archive=EntryArchive(m_context=context),
    )

    data.save(with_es=False)
    # assert archive
    assert test_archive.metadata.entry_type == 'PerovskiteSolarCell'
    assert test_archive.results.properties.optoelectronic.solar_cell.efficiency > -1
    assert len(test_archive.results.material.chemical_formula_reduced) > 0
    assert len(test_archive.data.jv.jv_curve[0].current_density) > -1
