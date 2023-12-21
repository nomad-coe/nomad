from math import factorial

# from posym import PointGroup, SymmetryObject
import numpy as np


class FillingRule:
    def __init__(self, occupation: int, orbital_indices: list[int]):
        self.occupation = occupation
        self.orbital_indices = orbital_indices

    def _fill(self, occupation: int, n_orbitals: int) -> list[bool]:
        raise NotImplementedError

    def fill(self) -> list[bool]:
        filling = self._fill(self.occupation, len(self.orbital_indices))
        occupations = np.zeros(32, dtype=bool)
        for index in sorted(self.orbital_indices):
            occupations[index] = filling.pop(0)
        return list(occupations)


class RussellSaundersState:
    @classmethod
    def generate_Js(cls, J1: float, J2: float, rising=True):
        J_min, J_max = sorted([abs(J1), abs(J2)])
        generator = range(
            int(J_max - J_min) + 1
        )  # works for both for fermions and bosons
        if rising:
            for jj in generator:
                yield J_min + jj
        else:
            for jj in generator:
                yield J_max - jj

    @classmethod
    def generate_MJs(cls, J, rising=True):
        generator = range(int(2 * J + 1))
        if rising:
            for m in generator:
                yield -J + m
        else:
            for m in generator:
                yield J - m

    def __init__(self, *args, **kwargs):
        self.J = kwargs.get('J')
        if self.J is None:
            raise TypeError
        self.occupation = kwargs.get('occ')
        if self.occupation is None:
            raise TypeError

    @property
    def multiplicity(self):
        return 2 * self.J + 1

    @property
    def degeneracy(self):
        return factorial(self.multiplicity) / (
            factorial(self.multiplicity - self.occupation) * factorial(self.occupation)
        )
