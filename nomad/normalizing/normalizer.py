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
from typing import List, Optional

from nomad.utils import get_logger
from nomad.metainfo import MSection
from nomad.datamodel import EntryArchive


class Normalizer(metaclass=ABCMeta):
    """
    A base class for normalizers. Only one instance of the Normalizer is
    created, and you should perform any heavy initialization in the constructor.
    The normalize method is called on archives, and this function should ideally
    not mutate the state of the shared normalizer instance.
    """

    domain: Optional[str] = 'dft'
    """Deprecated: The domain this normalizer should be used in. Default for all normalizer is 'DFT'."""
    normalizer_level = 0
    """Deprecated: Specifies the order of normalization with respect to other normalizers. Lower level
    is executed first."""

    def __init__(self, **kwargs) -> None:
        self.logger = get_logger(__name__)

    @abstractmethod
    def normalize(self, archive: EntryArchive, logger=None) -> None:
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)


class SystemBasedNormalizer(Normalizer, metaclass=ABCMeta):
    """
    A normalizer base class for normalizers that only touch a section_system.

    The normalizer is run on all section systems in a run. However, some systems,
    selected by heuristic, are more `representative systems` for the run. Sub-classes
    might opt to do additional work for the `representative systems`.

    Args:
        only_representatives: Will only normalize the `representative` systems.
    """

    def __init__(self, only_representatives: bool = False, **kwargs):
        super().__init__(**kwargs)
        self.only_representatives = only_representatives

    @property
    def quantities(self) -> List[str]:
        return [
            'atom_labels',
            'atom_positions',
            'atom_atom_numbers',
            'lattice_vectors',
            'simulation_cell',
            'configuration_periodic_dimensions',
        ]

    def _normalize_system(self, archive, system, is_representative):
        return self.normalize_system(archive, system, is_representative)

    @abstractmethod
    def normalize_system(
        self, archive: EntryArchive, system: MSection, is_representative: bool
    ) -> bool:
        """Normalize the given section and returns True, if successful"""
        pass

    def __representative_system(self, archive):
        """Used to select a representative system for this entry.

        Attempt to find a single section_system that is representative for the
        entry. The selection depends on the type of calculation.
        """
        system = None
        scc = None

        # Try to find workflow information and select the representative system
        # based on it
        workflow = archive.workflow2

        if workflow:
            try:
                iscc = workflow.results.calculation_result_ref
                system = iscc.system_ref
                if system is not None:
                    scc = iscc
            except Exception:
                pass

        # If no frame sequences detected, try to find valid scc by looping all
        # available in reverse order until a valid one is found.
        if system is None:
            try:
                sccs = archive.run[0].calculation
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
                    system = archive.run[0].system[-1]
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
            archive.run[0].m_cache['representative_scc_idx'] = scc.m_parent_index
        if system is not None:
            archive.run[0].m_cache['representative_system_idx'] = system.m_parent_index

        return system.m_resolved() if system is not None else None

    def __normalize_system(
        self, archive: EntryArchive, system, representative, logger=None
    ) -> bool:
        try:
            return self._normalize_system(archive, system, representative)

        except KeyError as e:
            self.logger.error(
                'could not read a system property',
                normalizer=self.__class__.__name__,
                section='system',
                g_index=system.m_parent_index,
                key_error=str(e),
                exc_info=e,
            )
            return False

        except Exception as e:
            self.logger.error(
                'Unexpected error during normalizing',
                normalizer=self.__class__.__name__,
                section='system',
                g_index=system.m_parent_index,
                exc_info=e,
                error=str(e),
            )
            raise e

    def normalize(self, archive: EntryArchive, logger=None) -> None:
        super().normalize(archive, logger)
        # If no section run detected, do nothing. Note that even the definition
        # may be missing and can cause an AttributeError.
        try:
            archive.run[0]
        except (AttributeError, IndexError):
            return

        # Process representative system first
        repr_sys_idx = None
        repr_sys = self.__representative_system(archive)
        if repr_sys is not None:
            repr_sys_idx = repr_sys.m_parent_index
            self.logger.info('chose "representative" section system')
            self.__normalize_system(archive, repr_sys, True, logger)

        # All the rest if requested
        if not self.only_representatives:
            for isys, system in enumerate(archive.run[0].system):
                if isys != repr_sys_idx:
                    self.__normalize_system(archive, system, False, logger)
