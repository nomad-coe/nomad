
class MObject:
    def __init__(self, m_definition: 'MElementDef', m_section: 'MSection' = None):
        self.m_definition = m_definition
        self.m_section = m_section

    @property
    def m_section(self):
        return self._section

    @m_section.setter
    def m_section(self, m_section: 'MSection'):
        self._section = m_section

        # add yourself to the parent section
        if m_section is not None:
            subsection = m_section.subsections.get(self.m_definition.name)
            if subsection is None:
                subsection = []
                m_section.subsections[self.m_definition.name] = subsection

            subsection.append(self)


class MSection(MObject):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.subsections = dict()
        self.properties = dict()

    def __getattr__(self, name: str):
        if name.startswith('new_'):
            name = name[len('new_'):]
            subsection = self.m_definition.get_subsection(name)

            if subsection is None:
                raise KeyError('Section "%s" does not have subsection "%s", available subsections are %s' % (self.m_definition.name, name, '?'))

            def constructor(**kwargs):
                new_section = subsection.impl(m_definition=subsection, m_section=self)
                for key, value in kwargs.items():
                    setattr(new_section, key, value)

                return new_section

            return constructor

    def __setattr__(self, name: str, value):
        try:
            super().__setattr__(name, value)
        except KeyError:
            self.__dict__[name] = value

    def __repr__(self):
        return ':%s' % self.m_definition.name


class MProperty(MObject):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.value = None


class MElementDef(MSection):
    def __init__(self, m_definition: 'MElementDef', name: str): # more **kwargs
        super().__init__(m_definition=m_definition)
        self.name = name

    def __repr__(self):
        return '%s:%s' % (self.name, self.m_definition.name)


class MSectionDef(MElementDef):
    def __init__(self, abstract: bool = False, repeats: bool = False, impl=MSection, **kwargs):
        super().__init__(**kwargs)
        self.repeats = repeats
        self.abstract = abstract
        self.impl = impl

    def get_subsection(self, name: str):
        return self.m_section.get_section(name)


m_section = MSectionDef(m_definition=None, name='section', repeats=True, impl=MSectionDef)
m_section.m_definition = m_section
m_element = MSectionDef(m_definition=m_section, name='element', abstract=True, impl=None)
m_package = MSectionDef(m_definition=m_section, name='package')


class MPackageDef(MElementDef):
    def __init__(self, name: str):
        super().__init__(m_definition=m_package, name=name)

    def get_section(self, name: str):
        # TODO add cache
        for section in self.subsections['section']:
            if section.name == name:
                return section

        return None


metainfo = MPackageDef(name='metainfo')
m_element.m_section = metainfo
m_section.m_section = metainfo
m_package.m_section = metainfo

m_property = metainfo.new_section(name='property', repeats=True, extends=m_element, section=m_package)
print(metainfo.subsections['section'])
metainfo.new_property(name='abstract', type=bool, section=m_element)
