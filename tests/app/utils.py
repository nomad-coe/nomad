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

from typing import List
import numpy as np

from nomad import processing, files
from nomad.datamodel import EntryMetadata, MongoMetadata, EntryArchive
from nomad.parsing import Backend

from tests.test_normalizing import run_normalize


class Upload():

    def __init__(self):
        self.entries: List[EntryMetadata] = []
        self.upload_id = 'test_upload_id'

    def create_upload_files(self) -> None:
        upload_files = files.StagingUploadFiles(self.upload_id, create=True)
        for entry_metadata in self.entries:
            archive = entry_metadata.m_parent
            if archive is None:
                archive = EntryArchive()
                archive.m_add_sub_section(EntryArchive.section_metadata, entry_metadata)

            upload_files.write_archive(entry_metadata.calc_id, archive.m_to_dict())

        upload_files.pack(self.entries, skip_raw=True)
        upload_files.delete()

        assert files.UploadFiles.get(self.upload_id) is not None

    def add_entry(self, entry_metadata: EntryMetadata):
        self.entries.append(entry_metadata)

        processing.Calc.create(
            calc_id=entry_metadata.calc_id,
            upload_id=entry_metadata.upload_id,
            mainfile=entry_metadata.mainfile,
            metadata=entry_metadata.m_to_dict(
                include_defaults=True, categories=[MongoMetadata])).save()

        entry_metadata.a_elastic.index()

    def create_test_structure(
            self, id: int, h: int, o: int, extra: List[str], periodicity: int,
            optimade: bool = True, metadata: dict = None):
        ''' Creates a calculation in Elastic and Mongodb with the given properties.

        Does require initialized :func:`elastic_infra` and :func:`mongo_infra`.

        Args:
            meta_info: A legace metainfo env.
            id: A number to create ``test_calc_id_<number>`` ids.
            h: The amount of H atoms
            o: The amount of O atoms
            extra: A list of further atoms
            periodicity: The number of dimensions to repeat the structure in
            optimade: A boolean. Iff true the entry will have optimade metadata. Default is True.
            metadata: Additional (user) metadata.
        '''
        atom_labels = ['H' for i in range(0, h)] + ['O' for i in range(0, o)] + extra
        test_vector = np.array([0, 0, 0])

        backend = Backend('public', False, True)  # type: ignore
        backend.openSection('section_run')
        backend.addValue('program_name', 'test_code')
        backend.openSection('section_system')

        backend.addArrayValues('atom_labels', np.array(atom_labels))
        backend.addArrayValues(
            'atom_positions', np.array([test_vector for i in range(0, len(atom_labels))]))
        backend.addArrayValues(
            'lattice_vectors', np.array([test_vector, test_vector, test_vector]))
        backend.addArrayValues(
            'configuration_periodic_dimensions',
            np.array([True for _ in range(0, periodicity)] + [False for _ in range(periodicity, 3)]))

        backend.closeSection('section_system', 0)
        backend.closeSection('section_run', 0)

        backend = run_normalize(backend)
        entry_metadata = backend.entry_archive.section_metadata

        entry_metadata.m_update(
            domain='dft', upload_id=self.upload_id, calc_id='test_calc_id_%d' % id,
            mainfile='test_mainfile', published=True, with_embargo=False)

        entry_metadata.apply_domain_metadata(backend)

        if metadata is not None:
            entry_metadata.m_update(**metadata)

        if not optimade:
            entry_metadata.dft.optimade = None

        self.add_entry(entry_metadata)
