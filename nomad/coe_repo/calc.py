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

from typing import List
import json
from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import relationship, aliased
from sqlalchemy.sql.expression import literal

from nomad import infrastructure, datamodel

from .user import User
from .base import Base, calc_citation_association, ownership, co_authorship, shareship, \
    Tag, Topics, CalcSet, calc_dataset_containment, Citation


class Calc(Base, datamodel.Calc):  # type: ignore
    __tablename__ = 'calculations'

    calc_id = Column(Integer, primary_key=True, autoincrement=True)
    origin_id = Column(Integer, ForeignKey('uploads.upload_id'))
    upload = relationship('Upload')
    checksum = Column(String)

    calc_meta_data = relationship('CalcMetaData', uselist=False, lazy='joined')
    user_meta_data = relationship('UserMetaData', uselist=False, lazy='joined')
    citations = relationship('Citation', secondary=calc_citation_association, lazy='joined')
    owners = relationship('User', secondary=ownership, lazy='joined')
    coauthors = relationship('User', secondary=co_authorship, lazy='joined')
    shared_with = relationship('User', secondary=shareship, lazy='joined')

    parents = relationship(
        'Calc',
        secondary=calc_dataset_containment,
        primaryjoin=calc_dataset_containment.c.children_calc_id == calc_id,
        secondaryjoin=calc_dataset_containment.c.parent_calc_id == calc_id,
        backref='children')

    @classmethod
    def load_from(cls, obj):
        repo_db = infrastructure.repository_db
        return repo_db.query(Calc).filter_by(calc_id=int(obj.pid)).first()

    @property
    def mainfile(self) -> str:
        return self.calc_meta_data.location

    @property
    def pid(self):
        return self.calc_id

    @property
    def comment(self) -> str:
        return self.user_meta_data.label

    @property
    def calc_hash(self) -> str:
        return self.checksum

    @property
    def references(self) -> List[str]:
        return list(citation.value for citation in self.citations if citation.kind == 'EXTERNAL')

    @property
    def uploader(self) -> User:
        assert len(self.owners) == 1, 'A calculation can only have one owner.'
        return self.owners[0]

    @property
    def with_embargo(self) -> bool:
        return self.user_meta_data.permission == 1

    @property
    def chemical_formula(self) -> str:
        return self.calc_meta_data.chemical_formula

    @property
    def filenames(self) -> List[str]:
        filenames = self.calc_meta_data.filenames.decode('utf-8')
        return json.loads(filenames)

    @property
    def all_datasets(self) -> List['DataSet']:
        assert self.calc_id is not None
        repo_db = infrastructure.repository_db
        query = repo_db.query(literal(self.calc_id).label('calc_id')).cte(recursive=True)
        right = aliased(query)
        left = aliased(CalcSet)
        query = query.union_all(repo_db.query(left.parent_calc_id).join(
            right, right.c.calc_id == left.children_calc_id))
        query = repo_db.query(query)
        dataset_calc_ids = list(r[0] for r in query if not r[0] == self.calc_id)
        if len(dataset_calc_ids) > 0:
            return [
                DataSet(dataset_calc)
                for dataset_calc in repo_db.query(Calc).filter(Calc.calc_id.in_(dataset_calc_ids))]
        else:
            return []

    @property
    def direct_datasets(self) -> List['DataSet']:
        return [DataSet(dataset_calc) for dataset_calc in self.parents]

    def set_value(self, topic_cid: int, value: str) -> None:
        if value is None:
            return

        repo_db = infrastructure.repository_db
        topic = repo_db.query(Topics).filter_by(topic=value).first()
        if not topic:
            topic = Topics(cid=topic_cid, topic=value)
            repo_db.add(topic)

        tag = Tag(calc=self, topic=topic)
        repo_db.add(tag)


class DataSet:
    def __init__(self, dataset_calc: Calc) -> None:
        self._dataset_calc = dataset_calc

    @property
    def id(self):
        return self._dataset_calc.calc_id

    @property
    def dois(self) -> List[Citation]:
        return list(citation for citation in self._dataset_calc.citations if citation.kind == 'INTERNAL')

    @property
    def name(self):
        return self._dataset_calc.calc_meta_data.chemical_formula
