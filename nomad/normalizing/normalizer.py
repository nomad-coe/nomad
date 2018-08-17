from abc import ABCMeta, abstractmethod

from nomad.parsing import AbstractParserBackend


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
