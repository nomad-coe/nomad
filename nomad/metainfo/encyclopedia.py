from elasticsearch_dsl import InnerDoc
from nomad.metainfo import MSection, Section, SubSection, Quantity


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


class Encyclopedia(MSection):
    m_def = Section(
        a_flask=dict(skip_none=True),
        a_elastic=dict(type=InnerDoc)
    )
    material = SubSection(sub_section=Material.m_def, repeats=False)
