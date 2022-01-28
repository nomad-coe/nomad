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

from abc import ABCMeta, abstractmethod
from typing import List

from nomad.utils import get_logger
from nomad.metainfo import MSection
from nomad.datamodel import EntryArchive


class Normalizer(metaclass=ABCMeta):
    '''
    A base class for normalizers. Normalizers work on a :class:`EntryArchive` section
    for read and write. Normalizer instances are reused.

    Arguments:
        entry_archive: The entry_archive root section of the archive to normalize.
    '''

    domain = 'dft'
    ''' The domain this normalizer should be used in. Default for all normalizer is 'DFT'. '''

    def __init__(self, entry_archive: EntryArchive) -> None:
        self.entry_archive = entry_archive
        try:
            self.section_run = entry_archive.run[0]
        except (AttributeError, IndexError):
            self.section_run = None
        self.logger = get_logger(__name__)

    @abstractmethod
    def normalize(self, logger=None) -> None:
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)


class SystemBasedNormalizer(Normalizer, metaclass=ABCMeta):
    '''
    A normalizer base class for normalizers that only touch a section_system.

    The normalizer is run on all section systems in a run. However, some systems,
    selected by heuristic, are more `representative systems` for the run. Sub-classes
    might opt to do additional work for the `representative systems`.

    Args:
        only_representatives: Will only normalize the `representative` systems.
    '''
    def __init__(self, entry_archive: EntryArchive, only_representatives: bool = False):
        super().__init__(entry_archive)
        self.only_representatives = only_representatives

    @property
    def quantities(self) -> List[str]:
        return [
            'atom_labels',
            'atom_positions',
            'atom_atom_numbers',
            'lattice_vectors',
            'simulation_cell',
            'configuration_periodic_dimensions'
        ]

    def _normalize_system(self, system, is_representative):
        return self.normalize_system(system, is_representative)

    @abstractmethod
    def normalize_system(self, system: MSection, is_representative: bool) -> bool:
        ''' Normalize the given section and returns True, iff successful'''
        pass

    def __representative_system(self):
        '''Used to select a representative system for this entry.

        Attempt to find a single section_system that is representative for the
        entry. The selection depends on the type of calculation.
        '''
        system = None
        scc = None

        # Try to find workflow information and select the representative system
        # based on it
        workflow = self.entry_archive.workflow

        if workflow:
            try:
                iscc = workflow.calculation_result_ref
                system = iscc.system_ref
                if system is not None:
                    scc = iscc
            except Exception:
                pass

        # If no frame sequences detected, try to find valid scc by looping all
        # available in reverse order until a valid one is found.
        if system is None:
            try:
                sccs = self.section_run.calculation
                for iscc in reversed(sccs):
                    isys = iscc.system_ref
                    if isys is not None:
                        system = isys
                        scc = iscc
                        break
            except Exception:
                pass

            # If no sccs exist, try to find systems
            if system is None:
                try:
                    system = self.section_run.system[-1]
                except Exception:
                    system = None

            if system is None:
                self.logger.warning('no "representative" section system found')

        # If the found system is referencing a subsystem, then we choose it as
        # the representative one. These smaller subsystems are much easier to
        # analyze. Currently used in phonon calculations.
        if system is not None:
            try:
                system_ref = system.sub_system_ref
                if system_ref is not None:
                    system = system_ref
            except Exception:
                pass

        if scc is not None:
            self.section_run.m_cache["representative_scc_idx"] = scc.m_parent_index
        if system is not None:
            self.section_run.m_cache["representative_system_idx"] = system.m_parent_index

        return system.m_resolved() if system is not None else None

    def __normalize_system(self, system, representative, logger=None) -> bool:
        try:
            return self._normalize_system(system, representative)

        except KeyError as e:
            self.logger.error(
                'could not read a system property', normalizer=self.__class__.__name__,
                section='system', g_index=system.m_parent_index, key_error=str(e), exc_info=e)
            return False

        except Exception as e:
            self.logger.error(
                'Unexpected error during normalizing', normalizer=self.__class__.__name__,
                section='system', g_index=system.m_parent_index, exc_info=e, error=str(e))
            raise e

    def normalize(self, logger=None) -> None:
        super().normalize(logger)
        # If no section run detected, do nothing
        if self.section_run is None:
            return

        # Process representative system first
        repr_sys_idx = None
        repr_sys = self.__representative_system()
        if repr_sys is not None:
            repr_sys_idx = repr_sys.m_parent_index
            self.logger.info('chose "representative" section system')
            self.__normalize_system(repr_sys, True, logger)

        # All the rest if requested
        if not self.only_representatives:
            for isys, system in enumerate(self.section_run.system):
                if isys != repr_sys_idx:
                    self.__normalize_system(system, False, logger)
