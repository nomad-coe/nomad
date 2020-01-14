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
s_sampling_method = "section_sampling_method"
r_scc_to_system = 'single_configuration_calculation_to_system_ref'
r_frame_sequence_local_frames = 'frame_sequence_local_frames_ref'
r_frame_sequence_to_sampling = "frame_sequence_to_sampling_ref"


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
        """Normalize the given section and returns True, if successful"""
        pass

    def __representative_systems(self):
        """Used to tag systems that are representative for a calculation. In
        practice, takes the first and two last frames of all
        section_frame_sequences. If no section_frame_sequence exists, take
        first and two last frames from all defined sccs.
        """
        systems = []
        sequences = []

        # Get all frame sequences
        try:
            frame_seqs = self._backend[s_frame_sequence]
        except Exception:
            frame_seqs = []
        else:
            for frame_seq in frame_seqs:
                try:
                    frames = frame_seq[r_frame_sequence_local_frames]
                except Exception:
                    pass
                else:
                    sequences.append(frames)

        # If no frames exist, consider all existing sccs
        if len(sequences) == 0:
            try:
                sccs = self._backend.get_sections(s_scc)
            except Exception:
                pass
            else:
                sequences.append(sccs)

        # Take first and last fwo frames from each detected sequence
        sccs = self._backend[s_scc]
        for seq in sequences:
            if len(seq) == 1:
                indices = [0]
            elif len(seq) == 2:
                indices = [0, -1]
            elif len(seq) > 2:
                indices = [0, -2, -1]
            else:
                break
            for scc_idx in [seq[idx] for idx in indices]:
                system_idx = sccs[scc_idx][r_scc_to_system]
                systems.append(system_idx)

        if len(systems) == 0:
            self.logger.error('no "representative" section system found')
        self.logger.info(
            'chose "representative" systems for normalization',
            number_of_systems=len(systems))

        return systems

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

        representative_systems = set(self.__representative_systems())
        all_systems = self._backend.get_sections(s_system)

        has_representative = False
        for g_index in representative_systems:
            has_representative = has_representative or self.__normalize_system(g_index, True, logger)

        # all the rest or until first representative depending on configuration
        if not self.only_representatives or not has_representative:
            for g_index in all_systems:
                if g_index not in representative_systems:
                    if self.__normalize_system(g_index, not has_representative, logger):
                        has_representative = True
                        if self.only_representatives:
                            break
