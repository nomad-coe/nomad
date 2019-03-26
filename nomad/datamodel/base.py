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

from typing import Iterable, List, Dict, Type
import datetime
from elasticsearch_dsl import Keyword

from nomad import utils, config


class UploadWithMetadata():
    """
    See :class:`CalcWithMetadata`.
    """

    def __init__(self, **kwargs):
        self.upload_id: str = None
        self.uploader: utils.POPO = None
        self.upload_time: datetime.datetime = None

        self.calcs: Iterable['CalcWithMetadata'] = list()

        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    def calcs_dict(self) -> Dict[str, 'CalcWithMetadata']:
        return {calc.calc_id: calc for calc in self.calcs}


class CalcWithMetadata():
    """
    A dict/POPO class that can be used for mapping calc representations with calc metadata.
    We have multi representations of calcs and their calc metadata. To avoid implement
    mappings between all combinations, just implement mappings with the class and use
    mapping transitivity. E.g. instead of A -> B, A -> this -> B.

    Attributes:
        upload_id: The ``upload_id`` of the calculations upload (random UUID).
        calc_id: The unique mainfile based calculation id.
        calc_hash: The raw file content based checksum/hash of this calculation.
        pid: The unique persistent id of this calculation.
        mainfile: The upload relative mainfile path.

        files: A list of all files, relative to upload.
        upload_time: The time when the calc was uploaded.
        uploader: An object describing the uploading user, has at least ``user_id``
        processed: Boolean indicating if this calc was successfully processed and archive
            data and calc metadata is available.
        last_processing: A datatime with the time of the last successful processing.
        nomad_version: A string that describes the version of the nomad software that was
            used to do the last successful processing.

        with_embargo: Show if user set an embargo on the calculation.
        coauthors: List of coauther user objects with at ``user_id``.
        shared_with: List of users this calcs ownership is shared with, objects with at ``user_id``.
        comment: String comment.
        references: Objects describing user provided references, keys are ``id`` and ``value``.
        datasets: Objects describing the datasets, keys are ``id``, ``name``, ``doi``.
            DOI is optional, is an object with key ``id``, ``value``.
    """
    def __init__(self, **kwargs):
        # id relevant metadata
        self.upload_id: str = None
        self.calc_id: str = None
        self.calc_hash: str = None
        self.mainfile: str = None
        self.pid: int = None

        # basic upload and processing related metadata
        self.upload_time: datetime.datetime = None
        self.files: List[str] = None
        self.uploader: utils.POPO = None
        self.processed: bool = False
        self.last_processing: datetime.datetime = None
        self.nomad_version: str = None

        # user metadata, i.e. quantities given and editable by the user
        self.with_embargo: bool = None
        self.published: bool = False
        self.coauthors: List[utils.POPO] = []
        self.shared_with: List[utils.POPO] = []
        self.comment: str = None
        self.references: List[utils.POPO] = []
        self.datasets: List[utils.POPO] = []

        # temporary reference to the backend after successful processing
        self.backend = None

        self.update(**kwargs)

    def to_dict(self):
        return {
            key: value for key, value in self.__dict__.items()
            if value is not None and key not in ['backend']
        }

    def update(self, **kwargs):
        for key, value in kwargs.items():
            if isinstance(value, list):
                if len(value) > 0 and isinstance(value[0], dict) and not isinstance(value[0], utils.POPO):
                    value = list(utils.POPO(**item) for item in value)
            if isinstance(value, dict) and not isinstance(value, utils.POPO):
                value = utils.POPO(**value)

            setattr(self, key, value)

    def apply_user_metadata(self, metadata: dict):
        """
        Applies a user provided metadata dict to this calc.
        """
        self.pid = metadata.get('_pid')
        self.comment = metadata.get('comment')
        self.upload_time = metadata.get('_upload_time')
        uploader_id = metadata.get('_uploader')
        if uploader_id is not None:
            self.uploader = utils.POPO(id=int(uploader_id))
        self.references = [utils.POPO(value=ref) for ref in metadata.get('references', [])]
        self.with_embargo = metadata.get('with_embargo', False)
        self.coauthors = [
            utils.POPO(id=int(user)) for user in metadata.get('coauthors', [])]
        self.shared_with = [
            utils.POPO(id=int(user)) for user in metadata.get('shared_with', [])]
        self.datasets = [
            utils.POPO(id=int(ds['id']), doi=utils.POPO(value=ds.get('_doi')), name=ds.get('_name'))
            for ds in metadata.get('datasets', [])]

    def apply_domain_metadata(self, backend):
        raise NotImplementedError()


class DomainQuantity:

    def __init__(
            self, description: str = None, multi: bool = False, aggregate: bool = False,
            elastic_mapping=None):

        self.name: str = None
        self.description = description
        self.multi = multi
        self.aggregate = aggregate
        self.elastic_mapping = elastic_mapping

        if self.elastic_mapping is None:
            self.elastic_mapping = Keyword(multi=self.multi)


class Domain:
    domain_class: Type[CalcWithMetadata] = None
    quantities: List[DomainQuantity] = []

    @classmethod
    def register_domain(cls, domain_class: type, domain_name: str, quantities: Dict[str, DomainQuantity]):
        assert cls.domain_class is None, 'you can only define one domain.'

        if not domain_name == config.domain:
            return

        cls.domain_class = domain_class

        reference_domain_calc = domain_class()
        reference_general_calc = CalcWithMetadata()

        for name, value in reference_domain_calc.__dict__.items():
            if not hasattr(reference_general_calc, name):
                quantity = quantities.get(name, None)
                if quantity is None:
                    quantity = DomainQuantity()
                    quantities[name] = quantity
                quantity.name = name
                quantity.multi = isinstance(value, list)
                cls.quantities.append(quantity)

        for name in quantities.keys():
            assert hasattr(reference_domain_calc, name) and not hasattr(reference_general_calc, name), \
                'quantity does not exist or overrides general non domain quantity'

        # utils.get_logger(__name__).info(
        #     'configured domain', domain=domain_name, domain_quantities=len(cls.quantities))
