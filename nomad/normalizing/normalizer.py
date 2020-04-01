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

from abc import ABCMeta, abstractmethod
from typing import List

from nomad.parsing import AbstractParserBackend
from nomad.utils import get_logger


class Normalizer(metaclass=ABCMeta):
    '''
    A base class for normalizers. Normalizers work on a :class:`AbstractParserBackend` instance
    for read and write. Normalizer instances are reused.

    Arguments:
        backend: The backend used to read and write data from and to.
    '''

    domain = 'dft'
    ''' The domain this normalizer should be used in. Default for all normalizer is 'DFT'. '''

    def __init__(self, backend: AbstractParserBackend) -> None:
        self._backend = backend
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
    def __init__(self, backend: AbstractParserBackend, only_representatives: bool = False):
        super().__init__(backend)
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

    def _normalize_system(self, g_index, is_representative):
        context = '/section_run/0/section_system/%d' % g_index

        self._backend.openContext(context)
        try:
            return self.normalize_system(g_index, is_representative)
        finally:
            self._backend.closeContext(context)

    @abstractmethod
    def normalize_system(self, section_system_index: int, is_representative: bool) -> bool:
        ''' Normalize the given section and returns True, iff successful'''
        pass

    def __representative_system(self):
        '''Used to select a representative system for this entry.

        Attempt to find a single section_system that is representative for the
        entry. The selection depends on the type of calculation.
        '''
        # Try to find a frame sequence, only first found is considered
        try:
            frame_seq = self._backend['section_frame_sequence'][-1]
            sec_sampling_method = frame_seq['frame_sequence_to_sampling_ref']
            sampling_method = sec_sampling_method['sampling_method']
            frames = frame_seq['frame_sequence_local_frames_ref']
            if sampling_method == 'molecular_dynamics':
                scc = frames[0]
            else:
                scc = frames[-1]
            return scc['single_configuration_calculation_to_system_ref']
        except Exception:
            pass

        # If no frame sequences detected, try use system of last scc
        try:
            sccs = self._backend['section_single_configuration_calculation']
            scc = sccs[-1]
            return scc['single_configuration_calculation_to_system_ref']
        except Exception:
            pass

        # If no sccs exist, use last system
        try:
            systems = self._backend['section_system']
            return systems[-1]
        except Exception:
            pass

        self.logger.error('no "representative" section system found')

        return None

    def __normalize_system(self, g_index, representative, logger=None) -> bool:
        try:
            return self._normalize_system(g_index, representative)

        except KeyError as e:
            self.logger.error(
                'could read a system property', normalizer=self.__class__.__name__,
                section='section_system', g_index=g_index, key_error=str(e), exc_info=e)
            return False

        except Exception as e:
            self.logger.error(
                'Unexpected error during normalizing', normalizer=self.__class__.__name__,
                section='section_system', g_index=g_index, exc_info=e, error=str(e))
            raise e

    def normalize(self, logger=None) -> None:
        super().normalize(logger)

        # Process representative system first
        representative_system = self.__representative_system()
        if representative_system is not None:
            representative_system_idx = representative_system.m_parent_index
            self.logger.info('chose "representative" section system')
            self.__normalize_system(representative_system_idx, True, logger)

        # All the rest if requested
        if not self.only_representatives:
            for system in self._backend['section_system']:
                if system != representative_system:
                    self.__normalize_system(system.m_parent_index, False, logger)
