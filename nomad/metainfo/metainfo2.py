from typing import Dict, Any


class MetainfoObjectMeta(type):

    def __new__(self, name, bases, dct):
        cls = super().__new__(self, name, bases, dct)

        for attr_name, attr in dct.items():
            attr_cls = attr.__class__
            if attr_cls.__module__ == __name__ and attr_cls.__name__ == 'Quantity':
                attr.init(cls, attr_name)
            if attr_cls.__module__ == __name__ and attr_cls.__name__ == 'Section':
                attr.init(cls)

        return cls


class MetainfoObject(metaclass=MetainfoObjectMeta):
    def __init__(self):
        self.m_data: Dict[str, Any] = {}


class BootstrapQuantity():
    def __init__(self, **kwargs):
        self.m_data = kwargs


class Quantity(MetainfoObject):
    # m_definition = Section(repeats=True, parent_section=Quantity)

    name = BootstrapQuantity(type=str)
    type = BootstrapQuantity(type=type)
    # section

    def __init__(self, **kwargs):
        super().__init__()
        name = kwargs.pop('name', None)
        if name is not None:
            self.m_data.update(name=name)

        for name, value in kwargs.items():
            setattr(self, name, value)

    def __get__(self, obj, type=None):
        return obj.m_data.get(self.m_data['name'])

    def __set__(self, obj, value):
        obj.m_data[self.m_data['name']] = value

    def init(self, cls, name):
        self.name = name


for name, attr in Quantity.__dict__.items():
    if isinstance(attr, BootstrapQuantity):
        quantity = Quantity(name=name, **attr.m_data)
        setattr(Quantity, name, quantity)
        quantity.init(Quantity, name)


class Section(MetainfoObject):
    # m_definition = Section(repeats=True, parent_section=Package)

    name = Quantity(type=str)
    repeats = Quantity(type=bool)
    # parent_section = Quantity(type=Reference(Section))
    # extends = Quantity(type=Reference(Section))
    # adds_to = Quantity(type=Reference(Section))

    def __init__(self, **kwargs):
        super().__init__()
        for name, value in kwargs.items():
            setattr(self, name, value)

    def __get__(self, obj, type=None):
        return self

    def init(self, cls):
        self.name = cls.__name__
        self.description = cls.__doc__.strip()


class System(MetainfoObject):
    """
    This is the system documentation
    """
    m_definition = Section(repeats=True)

    atom_label = Quantity(type=str, shape=['n'])


system = System()

system.atom_label = ['H', 'H', 'O']
print(system.m_data)
print(system.m_definition.m_data)

print(System.m_definition.name)
