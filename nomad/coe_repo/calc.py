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
from datetime import datetime

from nomad import infrastructure, utils
from nomad.datamodel import CalcWithMetadata

from . import base
from .user import User
from .base import Base, calc_citation_association, ownership, co_authorship, shareship, \
    Tag, Topics, CalcSet, calc_dataset_containment, Citation, Spacegroup, CalcMetaData, \
    CodeVersion, StructRatio, UserMetaData


class QueryCache:
    """
    Caches queries to avoid unnecessary flushes while bulk creating calcs.
    Faster than even SQLAlchemy with ``autoflush=False``, because of reasons.
    """

    def __init__(self):
        self._cache = {}

    def __call__(self, entity, **kwargs):
        key = json.dumps(dict(entity=entity.__class__.__name__, **kwargs))
        value = self._cache.get(key, None)
        if value is None:
            value = infrastructure.repository_db.query(entity).filter_by(**kwargs).first()
            if value is not None:
                self._cache[key] = value
        return value


class IllegalCalcMetadata(Exception): pass


class Calc(Base):
    __tablename__ = 'calculations'

    coe_calc_id = Column('calc_id', Integer, primary_key=True, autoincrement=True)
    origin_id = Column(Integer, ForeignKey('uploads.upload_id'))
    upload = relationship('Upload')
    checksum = Column(String)

    calc_metadata = relationship('CalcMetaData', uselist=False, lazy='joined')
    user_metadata = relationship('UserMetaData', uselist=False, lazy='joined')
    citations = relationship('Citation', secondary=calc_citation_association, lazy='joined')
    owners = relationship('User', secondary=ownership, lazy='joined')
    coauthors = relationship('User', secondary=co_authorship, lazy='joined')
    shared_with = relationship('User', secondary=shareship, lazy='joined')
    tags = relationship('Tag', lazy='subquery', join_depth=1)
    spacegroup = relationship('Spacegroup', lazy='joined', uselist=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.topic_ids = {}

    parents = relationship(
        'Calc',
        secondary=calc_dataset_containment,
        primaryjoin=calc_dataset_containment.c.children_calc_id == coe_calc_id,
        secondaryjoin=calc_dataset_containment.c.parent_calc_id == coe_calc_id,
        backref='children', lazy='subquery', join_depth=1)

    @classmethod
    def load_from(cls, obj):
        repo_db = infrastructure.repository_db
        return repo_db.query(Calc).filter_by(coe_calc_id=int(obj.pid)).first()

    @property
    def mainfile(self) -> str:
        return self.calc_metadata.location

    @property
    def pid(self) -> int:
        return self.coe_calc_id

    @property
    def comment(self) -> str:
        return self.user_metadata.label

    @property
    def calc_id(self) -> str:
        return self.checksum

    @property
    def references(self) -> List[str]:
        return list(citation.value for citation in self.citations if citation.kind == 'EXTERNAL')

    @property
    def uploader(self) -> User:
        assert len(self.owners) == 1, 'A calculation must have exactly one owner.'
        return self.owners[0]

    @property
    def with_embargo(self) -> bool:
        return self.user_metadata.permission == 1

    @property
    def formula(self) -> str:
        return self.calc_metadata.chemical_formula

    @property
    def files(self) -> List[str]:
        if self.calc_metadata is not None:
            if self.calc_metadata.filenames is not None:
                filenames = self.calc_metadata.filenames.decode('utf-8')
                return json.loads(filenames)

        return []

    @property
    def all_datasets(self) -> List['DataSet']:
        assert self.coe_calc_id is not None
        repo_db = infrastructure.repository_db
        query = repo_db.query(literal(self.coe_calc_id).label('coe_calc_id')).cte(recursive=True)
        right = aliased(query)
        left = aliased(CalcSet)
        query = query.union_all(repo_db.query(left.parent_calc_id).join(
            right, right.c.coe_calc_id == left.children_calc_id))
        query = repo_db.query(query)
        dataset_calc_ids = list(r[0] for r in query if not r[0] == self.coe_calc_id)
        if len(dataset_calc_ids) > 0:
            return [
                DataSet(dataset_calc)
                for dataset_calc in repo_db.query(Calc).filter(Calc.coe_calc_id.in_(dataset_calc_ids))]
        else:
            return []

    @property
    def direct_datasets(self) -> List['DataSet']:
        return [DataSet(dataset_calc) for dataset_calc in self.parents]

    def _set_value(self, topic_cid: int, value: str, cache: QueryCache) -> None:
        if value is None:
            return

        repo_db = infrastructure.repository_db
        topic = cache(Topics, cid=topic_cid, topic=value)
        if not topic:
            topic = Topics(cid=topic_cid, topic=value)
            repo_db.add(topic)
            repo_db.flush()

        if topic.tid not in self.topic_ids:
            tag = Tag(calc=self, topic=topic)
            self.topic_ids[topic.tid] = topic.tid
            repo_db.add(tag)
        else:
            logger = utils.get_logger(
                __name__, calc_id=self.calc_id, upload_id=self.upload.upload_id)
            logger.warning('double tag on same calc', cid=topic.cid, tid=topic.tid, value=topic.topic)

    _dataset_cache: dict = {}

    def apply_calc_with_metadata(self, calc: CalcWithMetadata, cache: QueryCache) -> None:
        """
        Applies the data from ``source`` to this coe Calc object.
        """
        repo_db = infrastructure.repository_db

        self.checksum = calc.calc_id
        source_code_version = calc.code_version  # TODO shorten version names
        code_version_obj = cache(CodeVersion, content=source_code_version)
        if code_version_obj is None:
            code_version_obj = CodeVersion(content=source_code_version)
            repo_db.add(code_version_obj)

        if calc.upload_time is not None:
            added_time = calc.upload_time
        elif self.upload is not None and self.upload.upload_time is not None:
            added_time = self.upload.upload_time
        else:
            added_time = datetime.now()

        metadata = CalcMetaData(
            calc=self,
            added=added_time,
            chemical_formula=calc.formula,
            filenames=('[%s]' % ','.join(['"%s"' % filename for filename in calc.files])).encode('utf-8'),
            location=calc.mainfile,
            version=code_version_obj)
        repo_db.add(metadata)

        struct_ratio = StructRatio(
            calc=self,
            chemical_formula=calc.formula,
            formula_units=1, nelem=len(calc.atoms))
        repo_db.add(struct_ratio)

        user_metadata = UserMetaData(
            calc=self,
            label=calc.comment,
            permission=(1 if calc.with_embargo else 0))
        repo_db.add(user_metadata)

        spacegroup = Spacegroup(calc=self, n=calc.spacegroup)
        repo_db.add(spacegroup)

        # topic based properties
        self._set_value(base.topic_code, calc.code_name, cache)
        for atom in set(calc.atoms):
            self._set_value(base.topic_atoms, str(atom), cache)
        self._set_value(base.topic_system_type, calc.system, cache)
        self._set_value(base.topic_xc_treatment, calc.xc_functional, cache)
        self._set_value(base.topic_crystal_system, calc.crystal_system, cache)
        self._set_value(base.topic_basis_set_type, calc.basis_set, cache)

        # user relations
        def add_users_to_relation(source_users, relation):
            for source_user in source_users:
                coe_user = cache(User, user_id=source_user.id)
                if coe_user is None:
                    raise IllegalCalcMetadata(
                        'User with user_id %s does not exist.' % source_user.id)
                source_user.update(coe_user.to_popo())
                relation.append(coe_user)

        if calc.uploader is not None:
            add_users_to_relation([calc.uploader], self.owners)
        elif self.upload is not None and self.upload.user is not None:
            self.owners.append(self.upload.user)
            calc.uploader = self.upload.user.to_popo()

        add_users_to_relation(calc.coauthors, self.coauthors)
        add_users_to_relation(calc.shared_with, self.shared_with)

        # datasets
        for dataset in calc.datasets:
            dataset_id = dataset.id
            coe_dataset_calc: Calc = cache(Calc, coe_calc_id=dataset_id)
            if coe_dataset_calc is None:
                coe_dataset_calc = Calc(coe_calc_id=dataset_id)
                repo_db.add(coe_dataset_calc)

                metadata = CalcMetaData(
                    calc=coe_dataset_calc,
                    added=self.upload.upload_time,
                    chemical_formula=dataset.name)
                repo_db.add(metadata)

                if dataset.doi is not None:
                    self._add_citation(coe_dataset_calc, dataset.doi['value'], 'INTERNAL', cache)

                # cause a flush to avoid future inconsistencies
                repo_db.flush()

            coe_dataset_rel = CalcSet(parent_calc_id=dataset_id, children_calc_id=self.coe_calc_id)
            repo_db.add(coe_dataset_rel)

            dataset.update(DataSet(coe_dataset_calc).to_popo())

        # references
        for reference in calc.references:
            self._add_citation(self, reference['value'], 'EXTERNAL', cache)

    def _add_citation(self, coe_calc: 'Calc', value: str, kind: str, cache: QueryCache) -> None:
        repo_db = infrastructure.repository_db
        citation = cache(Citation, value=value, kind=kind)

        if citation is None:
            citation = Citation(value=value, kind=kind)
            repo_db.add(citation)

        coe_calc.citations.append(citation)

    def to_calc_with_metadata(self) -> CalcWithMetadata:
        """
        Creates a :class:`CalcWithMetadata` instance with UCPM ids, and all UMD/CMD.
        Be aware that ``upload_id`` and ``calc_id``, might be old coe repository
        ``upload_name`` and calculation ``checksum`` depending on the context, i.e. used
        database.
        """
        result = CalcWithMetadata(
            upload_id=self.upload.upload_id if self.upload else None,
            calc_id=self.checksum)

        result.pid = self.pid
        result.mainfile = self.mainfile
        result.files = self.files

        for topic in [tag.topic for tag in self.tags]:
            if topic.cid == base.topic_code:
                result.code_name = topic.topic
            elif topic.cid == base.topic_basis_set_type:
                result.basis_set = topic.topic
            elif topic.cid == base.topic_xc_treatment:
                result.xc_functional = topic.topic
            elif topic.cid == base.topic_system_type:
                result.system = topic.topic
            elif topic.cid == base.topic_atoms:
                result.atoms.append(topic.topic)
            elif topic.cid == base.topic_crystal_system:
                result.crystal_system = topic.topic
            elif topic.cid in [1996, 1994, 703, 702, 701, 100]:
                # user/author, restriction, formulas?, another category
                pass
            else:
                raise KeyError('topic cid %s.' % str(topic.cid))

        result.code_version = self.calc_metadata.version.content
        result.formula = self.calc_metadata.chemical_formula
        result.spacegroup = self.spacegroup.n
        result.atoms.sort()

        datasets: List[DataSet] = []
        for parent in self.parents:
            parents = Calc._dataset_cache.get(parent, None)
            if parents is None:
                parents = parent.all_datasets
                Calc._dataset_cache[parent] = parents
            datasets.append(DataSet(parent))
            datasets.extend(parents)

        result.pid = self.pid
        result.uploader = self.uploader.to_popo()
        result.upload_time = self.calc_metadata.added
        result.datasets = list(ds.to_popo() for ds in datasets)
        result.with_embargo = self.with_embargo
        result.comment = self.comment
        result.references = list(
            citation.to_popo() for citation in self.citations
            if citation.kind == 'EXTERNAL')
        result.coauthors = list(user.to_popo() for user in self.coauthors)
        result.shared_with = list(user.to_popo() for user in self.shared_with)

        return result


class DataSet:
    def __init__(self, dataset_calc: Calc) -> None:
        self._dataset_calc = dataset_calc

    @property
    def id(self):
        return self._dataset_calc.coe_calc_id

    @property
    def doi(self) -> Citation:
        doi = None
        for citation in self._dataset_calc.citations:
            if citation.kind == 'INTERNAL':
                if doi is not None:
                    utils.get_logger(__name__).warning(
                        'dataset with multiple dois', dataset_id=self.id)
                doi = citation
        return doi

    @property
    def name(self):
        return self._dataset_calc.calc_metadata.chemical_formula

    def to_popo(self):
        return utils.POPO(id=self.id, doi=self.doi.to_popo(), name=self.name)
