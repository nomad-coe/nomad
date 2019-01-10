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
from sqlalchemy.orm import relationship

from nomad import infrastructure, datamodel

from .user import User
from .base import Base, calc_citation_association, ownership, co_authorship, shareship, \
    Tag, Topics


class Calc(Base, datamodel.Calc):  # type: ignore
    __tablename__ = 'calculations'

    calc_id = Column(Integer, primary_key=True, autoincrement=True)
    origin_id = Column(Integer, ForeignKey('uploads.upload_id'))
    upload = relationship('Upload')
    checksum = Column(String)

    calc_meta_data = relationship('CalcMetaData', uselist=False)
    user_meta_data = relationship('UserMetaData', uselist=False)
    citations = relationship('Citation', secondary=calc_citation_association)
    owners = relationship('User', secondary=ownership)
    coauthors = relationship('User', secondary=co_authorship)
    shared_with = relationship('User', secondary=shareship)

    @classmethod
    def create_from(cls, obj):
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
