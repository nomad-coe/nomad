from nomad.metainfo import MSection, Section, SubSection, Quantity
from nomad.metainfo.mongoengine_extension import Mongo, MongoDocument


class DOSSimilarity(MSection):
    m_def = Section(
        a_mongo=MongoDocument()
    )
    material_ids = Quantity(
        type=str,
        shape=["n_similar_materials"],
        a_mongo=Mongo(),
    )
    values = Quantity(
        type=float,
        shape=["n_similar_materials"],
        a_mongo=Mongo(),
    )


class Similarity(MSection):
    m_def = Section(
        a_mongo=MongoDocument()
    )
    electronic_dos = SubSection(sub_section=DOSSimilarity.m_def, repeats=False)


class Material(MSection):
    m_def = Section(
        a_mongo=MongoDocument()
    )
    material_id = Quantity(
        type=str,
        a_mongo=Mongo(primary_key=True)
    )
    similarity = SubSection(sub_section=Similarity.m_def, repeats=False)
