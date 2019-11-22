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
            self.normalize_system(g_index, is_representative)
        finally:
            self._backend.closeContext(context)

    @abstractmethod
    def normalize_system(self, section_system_index: int, is_representative: bool) -> None:
        """ Normalize the given section. """
        pass

    def __representative_systems(self):
        # look for sccs in last frames
        sccs = []
        try:
            frame_seqs = self._backend.get_sections(s_frame_sequence)
        except Exception:
            frame_seqs = []

        for frame_seq in frame_seqs:
            try:
                frames = self._backend.get_value(r_frame_sequence_local_frames, frame_seq)
            except Exception:
                frames = []

            if len(frames) > 0:
                sccs.append(frames[-1])

        # no sccs from frames -> consider all sccs
        if len(sccs) == 0:
            try:
                sccs = self._backend.get_sections(s_scc)
            except Exception:
                sccs = []

        try:
            systems = [self._backend.get_value(r_scc_to_system, scc) for scc in sccs]
        except Exception:
            systems = []

        # only take the first, and last two systems
        if len(systems) == 0:
            try:
                systems = self._backend.get_sections(s_system)
            except Exception:
                systems = []

        if len(systems) > 2:
            systems = [systems[0], systems[-2], systems[-1]]

        if len(systems) == 0:
            self.logger.error('no "representative" section system found')

        self.logger.info(
            'chose "representative" systems for normalization',
            number_of_systems=len(systems))

        return set(systems)

    def normalize(self, logger=None) -> None:
        super().normalize(logger)

        representative_systems = self.__representative_systems()
        all_systems = self._backend.get_sections(s_system)
        selected_systems = representative_systems if self.only_representatives else all_systems

        for g_index in selected_systems:
            try:
                self._normalize_system(g_index, g_index in representative_systems)
            except KeyError as e:
                self.logger.error(
                    'Could not read all input data', normalizer=self.__class__.__name__,
                    section='section_system', g_index=g_index, key_error=str(e))
            except Exception as e:
                self.logger.error(
                    'Unexpected error during normalizing', normalizer=self.__class__.__name__,
                    section='section_system', g_index=g_index, exc_info=e, error=str(e))
                raise e
