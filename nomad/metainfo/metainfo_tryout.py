"""
Some playground to try the CONCEPT.md ideas.
"""

from typing import Dict, List, Any, Union, Type
import json
import ase.data


class Units():

    def __getattribute__(self, name):
        return name


units = Units()


class Definition():
    m_definition: Any = None
    pass


class Property(Definition):
    pass


class Quantity(Property):
    def __init__(
            self,
            name: str = None,
            description: str = None,
            parent_section: 'Section' = None,
            shape: List[Union[str, int]] = [],
            type: Union['Enum', type] = None,
            unit: str = None,
            derived: bool = False,
            repeats: bool = False,
            synonym: str = None):

        self.name = name
        self.parent_section = parent_section.m_definition if parent_section is not None else None
        self.derived = derived
        self.synonym = synonym

    def __get__(self, obj: 'MetainfoObject', type: Type):
        if obj is None:
            # access on cls not obj
            return self

        if self.derived:
            derive_method = getattr(obj, 'm_derive_%s' % self.name, None)

            if derive_method is None:
                raise KeyError('Derived quantity %s is not implemented' % self.name)

            return derive_method()

        elif self.synonym is not None:
            return getattr(obj, self.synonym)

        else:
            return obj.m_data.get(self.name, None)

    def __set__(self, obj: 'MetainfoObject', value: Any):
        obj.m_data[self.name] = value

    def __delete__(self, obj: 'MetainfoObject'):
        del obj.m_data[self.name]

    def __repr__(self):
        base = self.name
        if self.parent_section is not None:
            return '%s.%s' % (str(self.parent_section), base)
        else:
            return base


class Section(Definition):
    def __init__(
            self,
            name: str = None,
            parent_section=None,
            repeats: bool = False,
            extends=None,
            adds_to=None):

        self.name = name
        self.parent_section = parent_section.m_definition if parent_section is not None else None
        self.repeats = repeats

    def __get__(self, obj: 'MetainfoObject', type: Type):
        return self

    def __repr__(self):
        base = self.name
        if self.parent_section is not None:
            return '%s.%s' % (str(self.parent_section), base)
        else:
            return base


class Reference(Property):
    pass


class Enum:
    def __init__(self, values: List[Any]):
        self.values = values


class MetainfoObjectMeta(type):
    def __new__(cls, cls_name, bases, dct):
        cls.m_definition = dct.get('m_definition', None)
        for name, value in dct.items():
            if isinstance(value, Property):
                value.name = name
                value.parent_section = cls.m_definition

        cls = super().__new__(cls, cls_name, bases, dct)

        if cls.m_definition is not None:
            if cls.m_definition.name is None:
                cls.m_definition.name = cls_name

        return cls


class MetainfoObject(metaclass=MetainfoObjectMeta):
    """
    Base class for all
    """
    m_definition: Any = None

    def __init__(self):
        self.m_data = dict(m_defintion=self.m_definition.name)

    def m_create(self, section_definition: Any, *args, **kwargs) -> Any:
        """
        Creates a sub section of the given section definition.
        """
        section_cls = section_definition
        definition = section_definition.m_definition
        sub_section = section_cls(*args, **kwargs)
        if definition.repeats:
            self.m_data.setdefault(definition.name, []).append(sub_section)
        else:
            # TODO test overwrite
            self.m_data[definition.name] = sub_section

        return sub_section

    def m_get_definition(self, name):
        """
        Returns the definition of the given quantity name.
        """
        descriptor = getattr(type(self), name)
        if descriptor is None:
            raise KeyError

        if not isinstance(descriptor, Property):
            raise KeyError

        return descriptor

    def m_to_dict(self) -> Dict[str, Any]:
        """
        Returns a JSON serializable dictionary version of this section (and all subsections).
        """
        return {
            key: value.m_to_dict() if isinstance(value, MetainfoObject) else value
            for key, value in self.m_data.items()
        }

    def m_to_json(self) -> str:
        return json.dumps(self.m_to_dict(), indent=2)

    def m_validate(self) -> bool:
        """
        Validates this sections content based on section and quantity definitions.
        Can be overwritten to customize the validation.
        """
        return True

    def __repr__(self) -> str:
        return self.m_to_json()


class Run(MetainfoObject):
    m_definition = Section()


class System(MetainfoObject):
    """
    The system is ...
    """
    m_definition = Section(parent_section=Run, repeats=True)

    n_atoms = Quantity(type=int, derived=True)

    atom_labels = Quantity(shape=['n_atoms'], type=Enum(ase.data.chemical_symbols))
    """
    Atom labels are ...
    """

    atom_species = Quantity(shape=['n_atoms'], type=int, derived=True)

    atom_positions = Quantity(shape=['n_atoms', 3], type=float, unit=units.m)

    cell = Quantity(shape=[3, 3], type=float, unit=units.m)
    lattice_vectors = Quantity(synonym='cell')

    pbc = Quantity(shape=[3], type=bool)

    def m_derive_atom_species(self) -> List[int]:
        return [ase.data.atomic_numbers[label] for label in self.atom_labels]

    def m_derive_n_atoms(self) -> int:
        return len(self.atom_labels)


class VaspSystem(MetainfoObject):
    m_definition = Section(adds_to=System)

    vasp_specific = Quantity(type=str)


run = Run()
print(run.m_definition)

system = run.m_create(System)
system.atom_labels = ['H', 'H', 'O']
system.atom_positions = [[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0]]
system.cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
system.pbc = [False, False, False]

print(system.atom_species)
print(system.lattice_vectors)
print(system.n_atoms)

print(system.__class__.m_definition)
print(system.m_definition)
print(system.m_get_definition('atom_labels'))

print(run)
