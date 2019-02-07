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

"""
Declarative SQLAlchemy base model definitions for the repository db schema. Does
not include the *big* datamodel entities: `User`, `Upload`, `Calc`; they can
be found in their own submodules.
"""

from sqlalchemy import Column, Integer, String, DateTime, ForeignKey, Enum, Table
from sqlalchemy.orm import relationship
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import BYTEA
from sqlalchemy.ext.declarative import declarative_base


Base = declarative_base()

calc_citation_association = Table(
    'metadata_citations', Base.metadata,
    Column('calc_id', Integer, ForeignKey('calculations.calc_id')),
    Column('citation_id', Integer, ForeignKey('citations.citation_id')))


ownership = Table(
    'ownerships', Base.metadata,
    Column('calc_id', Integer, ForeignKey('calculations.calc_id')),
    Column('user_id', Integer, ForeignKey('users.user_id')))

co_authorship = Table(
    'coauthorships', Base.metadata,
    Column('calc_id', Integer, ForeignKey('calculations.calc_id')),
    Column('user_id', Integer, ForeignKey('users.user_id')))

shareship = Table(
    'shareships', Base.metadata,
    Column('calc_id', Integer, ForeignKey('calculations.calc_id')),
    Column('user_id', Integer, ForeignKey('users.user_id')))


class CalcMetaData(Base):  # type: ignore
    __tablename__ = 'metadata'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    added = Column(DateTime)
    chemical_formula = Column(String)
    filenames = Column(BYTEA)
    location = Column(String)
    version_id = Column(Integer, ForeignKey('codeversions.version_id'))
    version = relationship('CodeVersion', lazy='joined', uselist=False)


class UserMetaData(Base):  # type: ignore
    __tablename__ = 'user_metadata'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    label = Column(String)
    calc = relationship('Calc')
    permission = Column(Integer)


class StructRatio(Base):  # type: ignore
    __tablename__ = 'struct_ratios'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    formula_units = Column(Integer)
    nelem = Column(Integer)
    chemical_formula = Column(String)


class CodeVersion(Base):  # type: ignore
    __tablename__ = 'codeversions'

    version_id = Column(Integer, primary_key=True, autoincrement=True)
    content = Column(String)


class Spacegroup(Base):  # type: ignore
    __tablename__ = 'spacegroups'

    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    n = Column(Integer)


topic_code = 220
topic_atoms = 10
topic_system_type = 50
topic_xc_treatment = 75
topic_crystal_system = 90
topic_basis_set_type = 80


class Tag(Base):  # type: ignore
    __tablename__ = 'tags'
    calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    calc = relationship('Calc')
    tid = Column(Integer, ForeignKey('topics.tid'), primary_key=True)
    topic = relationship('Topics', lazy='joined', uselist=False)

    def __repr__(self):
        return '<Tag(calc_id="%d", tid="%d)>' % (int(self.calc_id), int(self.tid))


class Topics(Base):  # type: ignore
    __tablename__ = 'topics'
    tid = Column(Integer, primary_key=True, autoincrement=True)
    cid = Column(Integer)
    topic = Column(String)


class CalcSet(Base):  # type: ignore
    __tablename__ = 'calcsets'

    parent_calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)
    children_calc_id = Column(Integer, ForeignKey('calculations.calc_id'), primary_key=True)


calc_dataset_containment = Table(
    'calcsets', Base.metadata, extend_existing=True)


class Citation(Base):  # type: ignore
    __tablename__ = 'citations'

    citation_id = Column(Integer, primary_key=True)
    value = Column(String)
    kind = Column(Enum('INTERNAL', 'EXTERNAL', name='citation_kind_enum'))

    def to_dict(self) -> dict:
        return dict(id=self.citation_id, value=self.value)
