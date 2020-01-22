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

s_system = 'section_system'
s_scc = 'section_single_configuration_calculation'
s_frame_sequence = 'section_frame_sequence'
r_scc_to_system = 'single_configuration_calculation_to_system_ref'
r_frame_sequence_local_frames = 'frame_sequence_local_frames_ref'


class Normalizer(metaclass=ABCMeta):
    """
    A base class for normalizers. Normalizers work on a :class:`AbstractParserBackend` instance
    for read and write. Normalizer instances are reused.

    Arguments:
        backend: The backend used to read and write data from and to.
    """

    domain = 'DFT'
    """ The domain this normalizer should be used in. Default for all normalizer is 'DFT'. """

    def __init__(self, backend: AbstractParserBackend) -> None:
        self._backend = backend
        self.logger = get_logger(__name__)

    @abstractmethod
    def normalize(self, logger=None) -> None:
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
        """ Normalize the given section and returns True, iff successful"""
        pass

    def __representative_system(self):
        """Used to select a representative system for this entry."""

        system_idx = None

        # Try to find a frame sequence, only first found is considered
        try:
            frame_seqs = self._backend[s_frame_sequence]
            frame_seq = frame_seqs[0]
            sampling_method_idx = frame_seq["frame_sequence_to_sampling_ref"]
            sec_sampling_method = self._backend["section_sampling_method"][sampling_method_idx]
            sampling_method = sec_sampling_method["sampling_method"]
            frames = frame_seq["frame_sequence_local_frames_ref"]
            if sampling_method == "molecular_dynamics":
                scc_idx = frames[0]
            else:
                scc_idx = frames[-1]
            scc = self._backend[s_scc][scc_idx]
            system_idx = scc["single_configuration_calculation_to_system_ref"]
        except Exception:
            frame_seqs = []

        # If no frame sequences detected, try to find scc
        if len(frame_seqs) == 0:
            try:
                sccs = self._backend[s_scc]
                scc = sccs[-1]
                system_idx = scc["single_configuration_calculation_to_system_ref"]
            except KeyError:
                sccs = []

            # If no sccs exist, try to find systems
            if len(sccs) == 0:
                try:
                    systems = self._backend.get_sections(s_system)
                    system = systems[-1]
                    system_idx = system["single_configuration_calculation_to_system_ref"]
                except KeyError:
                    sccs = []

            if system_idx is None:
                self.logger.error('no "representative" section system found')
            else:
                self.logger.info(
                    'chose "representative" system for normalization',
                )

        return system_idx

    def __normalize_system(self, g_index, representative, logger=None) -> bool:
        try:
            return self._normalize_system(g_index, representative)

        except KeyError as e:
            self.logger.error(
                'Could not read all input data', normalizer=self.__class__.__name__,
                section='section_system', g_index=g_index, key_error=str(e))
            return False

        except Exception as e:
            self.logger.error(
                'Unexpected error during normalizing', normalizer=self.__class__.__name__,
                section='section_system', g_index=g_index, exc_info=e, error=str(e))
            raise e

    def normalize(self, logger=None) -> None:
        super().normalize(logger)

        all_systems = self._backend.get_sections(s_system)

        # Process representative system first
        representative_system_idx = self.__representative_system()
        if representative_system_idx is not None:
            self.__normalize_system(representative_system_idx, True, logger)

        # All the rest if requested
        if not self.only_representatives:
            for g_index in all_systems:
                if g_index != representative_system_idx:
                    self.__normalize_system(g_index, False, logger)
