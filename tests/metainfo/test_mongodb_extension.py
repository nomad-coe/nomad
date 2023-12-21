#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import json
import numpy as np
from nomad.metainfo import MSection, Section, Quantity, SubSection
from nomad.metainfo.mongoengine_extension import MongoDocument, Mongo


class B(MSection):
    m_def = Section(a_mongo=MongoDocument())
    value = Quantity(type=str, a_mongo=Mongo())


class C(MSection):
    m_def = Section(a_mongo=MongoDocument())
    value = Quantity(type=str, a_mongo=Mongo())


class D(MSection):
    m_def = Section()
    value = Quantity(type=str, a_mongo=Mongo())


class A(MSection):
    """Root level document with primary key."""

    m_def = Section(a_mongo=MongoDocument())
    primary_id = Quantity(type=str, a_mongo=Mongo(primary_key=True))
    array = Quantity(type=float, shape=[2], a_mongo=Mongo())
    not_in_mongo = Quantity(type=str)
    value1 = Quantity(type=int, a_mongo=Mongo())
    value2 = Quantity(type=int, a_mongo=Mongo())
    b = SubSection(sub_section=B.m_def)
    c = SubSection(sub_section=C.m_def, repeats=True)
    d = SubSection(sub_section=D.m_def)


def test_create_new(mongo):
    a = A()
    a.primary_id = '123'
    a.not_in_mongo = 'not_in_mongo'
    b = a.m_create(B)
    b.value = 'b_value'
    c = a.m_create(C)
    c.value = 'c_value'
    c = a.m_create(C)
    c.value = 'c_value'

    # Create JSON with the values that are supposed to be in mongo
    a_dict = a.m_to_dict()
    del a_dict['not_in_mongo']
    a_json = json.dumps(a_dict, sort_keys=True)

    # Store to mongo
    mongo_doc = a.a_mongo
    mongo_doc.save()

    # Retrieve from mongo, and convert to JSON
    a_from_db = A.m_def.a_mongo.get(primary_id='123')
    a_from_db_json = json.dumps(a_from_db.m_to_dict(), sort_keys=True)

    # Test equality of the JSON serializations
    assert a_json == a_from_db_json


def test_update_with_new():
    a = A()
    a.primary_id = '123'
    a.value1 = 1
    a.value2 = 2

    # Store to mongo
    a.a_mongo.save()

    # Update with new document that has the same ID
    a_new = A()
    a_new.primary_id = '123'
    a_new.value2 = 3
    a_new.a_mongo.save()

    # Check that the document has only partly been updated
    a_from_db = A.m_def.a_mongo.get(primary_id='123')
    assert a_from_db.value1 == 1
    assert a_from_db.value2 == 3


def test_update_self():
    a = A()
    a.primary_id = '123'
    a.value1 = 1
    a.value2 = 2

    # Store to mongo
    a.a_mongo.save()

    # Update the metainfo and resave
    a.value2 = 3
    a.a_mongo.save()

    # Check that the document has only partly been updated
    a_from_db = A.m_def.a_mongo.get(primary_id='123')
    assert a_from_db.value1 == 1
    assert a_from_db.value2 == 3


def test_annotations(mongo):
    """Test that non-annotated quantities and sections are not stored."""
    a = A()
    a.primary_id = '123'
    a.not_in_mongo = 'not_in_mongo'
    d = a.m_create(D)
    d.value = 'b_value'

    # Store to mongo
    a.a_mongo.save()

    # Check that values do not exist in mongodb
    a_from_db = A.m_def.a_mongo.get(primary_id='123')
    assert a_from_db.not_in_mongo is None
    assert a_from_db.d is None


def test_repeated_subsections():
    a = A()
    a.primary_id = '123'
    c = a.m_create(C)
    c.value = 'c_value'
    c = a.m_create(C)
    c.value = 'c_value'

    # Store to mongo
    a.a_mongo.save()

    # Check that both sections are stored in mongodb
    a_from_db = A.m_def.a_mongo.get(primary_id='123')
    assert len(a_from_db.c) == 2


def test_arrays():
    a = A()
    a.primary_id = '123'
    a.array = np.array([1.2, 3.4])

    # Store to mongo
    a.a_mongo.save()

    # Check that array was correctly save
    a_from_db = A.m_def.a_mongo.get(primary_id='123')
    assert a_from_db.array == [1.2, 3.4]
