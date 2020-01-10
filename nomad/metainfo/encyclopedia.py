from elasticsearch_dsl import InnerDoc
from nomad.metainfo import MSection, Section, SubSection, Quantity, Enum


class Material(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Section for storing the data that links this entry into a specific material.
        """
    )
    material_hash = Quantity(
        type=str,
        description="""
        A unique material identifier. For crystals the hash
        identifier is constructed from formula, space group and
        wyckoff_position_population.
        """
    )


class Calculation(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc),
        description="""
        Section for storing data related to a calculation that is identified
        from this entry.
        """
    )
    run_type = Quantity(
        type=Enum("single point", "geometry optimization", "molecular dynamics", "phonon calculation", "QHA calculation", "GW calculation", "equation of state", "parameter variation", "unavailable"),
        description="""
        Defines the type of run identified for this entry.
        """
    )


class Encyclopedia(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc)
    )
    material = SubSection(sub_section=Material.m_def, repeats=False)
    calculation = SubSection(sub_section=Calculation.m_def, repeats=False)
