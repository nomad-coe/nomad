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

from typing import List, Dict, Any
import json
from sqlalchemy import Column, Integer, String, ForeignKey
from sqlalchemy.orm import relationship, aliased
from sqlalchemy.sql.expression import literal
from datetime import datetime
import os.path

from nomad import infrastructure, utils, config, files
from nomad.datamodel import DFTCalcWithMetadata

from . import base
from .user import User
from .base import Base, calc_citation_association, ownership, co_authorship, shareship, \
    Tag, Topics, CalcSet, calc_dataset_containment, Citation, Spacegroup, CalcMetaData, \
    CodeVersion, StructRatio, UserMetaData


handle_base = '0123456789abcdefghijklmnopqrstuvwxyz'


def create_handle(pid: int) -> str:
    """
    Create a handle for the given pid. The pid is an autoincrement number. The handle
    a 'base32' encoded string of that number. Therefore, its string representation is a
    little shorter. The handle is prefixed with the configured handle prefix.
    """

    value = pid
    result = ''
    while value > 0:
        result += handle_base[value & 31]
        value = value >> 5

    return config.repository_db.handle_prefix + result[::-1]


class PublishContext:
    """
    Utilities necessary during adding calculations to the repo db.
    Caches queries to avoid unnecessary flushes while bulk creating calcs.
    Faster than even SQLAlchemy with ``autoflush=False``, because of reasons.
    Access to a logger with bound data about the upload, etc.
    """

    def __init__(self, upload_id: str = None, **kwargs):
        self._cache: Dict[str, Any] = {}
        self.upload_id = upload_id
        self.upload_files = None if upload_id is None else files.UploadFiles.get(upload_id, is_authorized=lambda: True)
        self.logger = utils.get_logger(__name__, upload_id=upload_id, **kwargs)

    def cache(self, entity, **kwargs):
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
    handlepid = Column(String)
    origin_id = Column(Integer, ForeignKey('uploads.upload_id'))
    upload = relationship('Upload', lazy='joined')
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

    @staticmethod
    def from_calc_id(calc_id: str) -> 'Calc':
        repo_db = infrastructure.repository_db
        calcs = repo_db.query(Calc).filter_by(checksum=calc_id)
        assert calcs.count() <= 1, 'Calc id/checksum must be unique'
        return calcs.first()

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
        # permission = 1 means public
        # permission = 0 means not public, i.e. with embargo
        return self.user_metadata.permission != 1

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

    def _set_value(self, topic_cid: int, value: str, context: PublishContext) -> None:
        if value is None:
            return

        repo_db = infrastructure.repository_db
        topic = context.cache(Topics, cid=topic_cid, topic=value)
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

    def apply_calc_with_metadata(self, calc: DFTCalcWithMetadata, context: PublishContext) -> None:
        """
        Applies the data from ``source`` to this coe Calc object.
        """
        repo_db = infrastructure.repository_db

        self.checksum = calc.calc_id
        source_code_version = calc.code_version  # TODO shorten version names
        code_version_obj = context.cache(CodeVersion, content=source_code_version)
        if code_version_obj is None:
            code_version_obj = CodeVersion(content=source_code_version)
            repo_db.add(code_version_obj)
            repo_db.flush()

        if calc.upload_time is not None:
            added_time = calc.upload_time
        elif self.upload is not None and self.upload.upload_time is not None:
            added_time = self.upload.upload_time
        else:
            added_time = datetime.utcnow()

        upload_id = context.upload_id
        upload_files = context.upload_files
        coe_files = list()
        if upload_files is None:
            upload_size = -1
        else:
            upload_size = 0

        for calc_file in calc.files:
            if config.repository_db.mode == 'coe':
                coe_file = os.path.join('$EXTRACTED', 'nomad', upload_id, calc_file).replace('/', '\\/')
            else:
                coe_file = calc_file

            if upload_files is not None:
                upload_size += upload_files.raw_file_size(calc_file)
            coe_files.append(coe_file)

        metadata = CalcMetaData(
            calc=self,
            added=added_time,
            oadate=added_time,
            chemical_formula=calc.formula,
            filenames=('[%s]' % ','.join(['"%s"' % coe_file for coe_file in coe_files])).encode('utf-8'),
            download_size=upload_size,
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
            permission=(0 if calc.with_embargo else 1))
        repo_db.add(user_metadata)

        if isinstance(calc.spacegroup, int) or calc.spacegroup.isdigit():
            spacegroup = Spacegroup(calc=self, n=calc.spacegroup)
        else:
            spacegroup = Spacegroup(calc=self, n='0')
        repo_db.add(spacegroup)

        # topic based properties
        self._set_value(base.topic_code, calc.code_name, context)
        for atom in set(calc.atoms):
            self._set_value(base.topic_atoms, str(atom), context)
        self._set_value(base.topic_system_type, calc.system, context)
        self._set_value(base.topic_xc_treatment, calc.xc_functional, context)
        self._set_value(base.topic_crystal_system, calc.crystal_system, context)
        self._set_value(base.topic_basis_set_type, calc.basis_set, context)

        # user relations
        def add_users_to_relation(source_users, relation):
            for source_user in source_users:
                coe_user = context.cache(User, user_id=int(source_user.id))
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
        calcs_existing_datasets: List[int] = []
        for dataset in calc.datasets:
            dataset_id = dataset.id
            if dataset_id in calcs_existing_datasets:
                continue
            else:
                calcs_existing_datasets.append(dataset_id)

            coe_dataset_calc: Calc = context.cache(Calc, coe_calc_id=dataset_id)
            if coe_dataset_calc is None:
                coe_dataset_calc = Calc(coe_calc_id=dataset_id)
                repo_db.add(coe_dataset_calc)

                metadata = CalcMetaData(
                    calc=coe_dataset_calc,
                    added=self.upload.upload_time,
                    chemical_formula=dataset.name)
                repo_db.add(metadata)
                repo_db.flush()

                if dataset.doi is not None:
                    self._add_citation(coe_dataset_calc, dataset.doi['value'], 'INTERNAL', context)

                # cause a flush to create the backdirection of the above established
                # metadata-dataset_calc relation
                repo_db.flush()

            self.parents.append(coe_dataset_calc)

            dataset.update(DataSet(coe_dataset_calc).to_popo())

        # references
        for reference in calc.references:
            self._add_citation(self, reference['value'], 'EXTERNAL', context)

        repo_db.flush()

    def _add_citation(self, coe_calc: 'Calc', value: str, kind: str, context: PublishContext) -> None:
        if value is None or kind is None:
            context.logger.warning(
                'citation without value or kind str', value=value, kind=kind, calc_id=self.calc_id)
            return

        repo_db = infrastructure.repository_db
        citation = context.cache(Citation, value=value, kind=kind)

        if citation is None:
            citation = Citation(value=value, kind=kind)
            repo_db.add(citation)

        coe_calc.citations.append(citation)

    def to_calc_with_metadata(self) -> DFTCalcWithMetadata:
        """
        Creates a :class:`DFTCalcWithMetadata` instance with UCPM ids, and all UMD/CMD.
        Be aware that ``upload_id`` and ``calc_id``, might be old coe repository
        ``upload_name`` and calculation ``checksum`` depending on the context, i.e. used
        database.
        """
        result = DFTCalcWithMetadata(
            upload_id=self.upload.upload_id if self.upload else None,
            calc_id=self.checksum)

        result.pid = self.pid
        result.mainfile = self.mainfile
        result.files = self.files

        for topic in [tag.topic for tag in self.tags]:
            if topic is None:
                continue

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
        if self.spacegroup is not None:
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
        return utils.POPO(
            id=self.id,
            name=self.name,
            doi=self.doi.to_popo() if self.doi is not None else None)
