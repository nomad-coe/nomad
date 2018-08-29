from abc import ABCMeta, abstractmethod
from typing import List, Dict, Any

from nomad.parsing import AbstractParserBackend
from nomad.utils import get_logger

logger = get_logger(__name__)


class Normalizer(metaclass=ABCMeta):
    """
    A base class for normalizers. Normalizers work on a :class:`AbstractParserBackend` instance
    for read and write.
    """
    def __init__(self, backend: AbstractParserBackend) -> None:
        self._backend = backend

    @abstractmethod
    def normalize(self) -> None:
        pass


class SystemBasedNormalizer(Normalizer, metaclass=ABCMeta):

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

    def _normalize_system(self, g_index):
        input_data = dict(
            uri='/section_run/0/section_system/%d' % g_index,
            gIndex=g_index)
        for quantity in self.quantities:
            try:
                input_data[quantity] = self._backend.get_value(quantity, g_index)
            except KeyError:
                # onyl fail when the normalizer actually uses the respecitive value
                pass

        context = input_data['uri']
        self._backend.openContext(context)
        try:
            self.normalize_system(input_data)
        finally:
            self._backend.closeContext(context)

    @abstractmethod
    def normalize_system(self, section_system: Dict[str, Any]) -> None:
        pass

    def normalize(self) -> None:
        for g_index in self._backend.get_sections('section_system'):
            try:
                self._normalize_system(g_index)
            except KeyError as e:
                logger.error(
                    'Could not read all input data', normalizer=self.__class__.__name__,
                    section='section_system', g_index=g_index, key_error=str(e))
            except Exception as e:
                logger.error(
                    'Unexpected error during normalizing', normalizer=self.__class__.__name__,
                    section='section_system', g_index=g_index, exc_info=e)
                raise e
