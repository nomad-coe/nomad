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

from elasticsearch_dsl import Document

from nomad.metainfo import MSection, Section, Quantity, SubSection, Reference
from nomad.metainfo.elastic_extension import ElasticDocument, Elastic


class User(MSection):

    user_id = Quantity(type=str, a_elastic=Elastic())
    name = Quantity(type=str, a_elastic=Elastic())
    email = Quantity(type=str)


class EMS(MSection):
    experiment_method = Quantity(type=str)


class DFT(MSection):
    atoms = Quantity(
        type=str, shape=['0..*'],
        a_elastic=[Elastic(), Elastic()])
    n_atoms = Quantity(
        type=int, derived=lambda x: len(x.atoms),
        a_elastic=Elastic(index=False, field='natoms'))
    only_atoms = Quantity(
        type=str, shape=['0..*'], derived=lambda x: x.atoms,
        a_elastic=Elastic(value=lambda x: ','.join(x.atoms)))


class Entry(MSection):
    m_def = Section(a_elastic=ElasticDocument(id=lambda entry: entry.entry_id))

    entry_id = Quantity(type=str, a_elastic=Elastic())
    uploader = Quantity(type=Reference(User.m_def), a_elastic=Elastic())

    dft = SubSection(sub_section=DFT.m_def)
    ems = SubSection(sub_section=EMS.m_def)


def test_document():
    annotation = Entry.m_def.a_elastic
    document = annotation.document

    assert issubclass(document, Document)
    assert DFT.m_def.qualified_name() in ElasticDocument._all_documents
    assert EMS.m_def.qualified_name() not in ElasticDocument._all_documents

    assert DFT.atoms.a_elastic[0].qualified_field == 'dft.atoms'
    assert DFT.n_atoms.a_elastic.qualified_field == 'dft.natoms'


def test_create_entry():
    user = User(user_id='test_user', name='Test Tester', email='test@testing.com')
    entry = Entry(entry_id='test_id', uploader=user)
    entry.m_create(DFT).atoms = ['H', 'O']

    index_entry = Entry.m_def.a_elastic.create_index_entry(entry)

    assert index_entry.entry_id == entry.entry_id
    assert index_entry.uploader.name == 'Test Tester'
    assert index_entry.dft.atoms == ['H', 'O']
    assert index_entry.dft.natoms == 2
    assert index_entry.dft.only_atoms == 'H,O'
    assert hasattr(index_entry, 'ems') is False
