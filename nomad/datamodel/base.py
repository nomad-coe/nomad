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

from typing import Iterable, List, Dict, Type, Tuple, Callable, Any
import datetime
from elasticsearch_dsl import Keyword
from collections.abc import Mapping
import numpy as np

from nomad import utils, config
from nomad.metainfo import MSection


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


class CalcWithMetadata(Mapping):
    """
    A dict/POPO class that can be used for mapping calc representations with calc metadata.
    We have multi representations of calcs and their calc metadata. To avoid implement
    mappings between all combinations, just implement mappings with the class and use
    mapping transitivity. E.g. instead of A -> B, A -> this -> B.

    This is basically an abstract class and it has to be subclassed for each :class:`Domain`.
    Subclasses can define additional attributes and have to implement :func:`apply_domain_metadata`
    to fill these attributes from processed entries, i.e. instance of :class:`nomad.parsing.LocalBackend`.

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
        self.raw_id: str = None

        # basic upload and processing related metadata
        self.upload_time: datetime.datetime = None
        self.files: List[str] = None
        self.uploader: utils.POPO = None
        self.processed: bool = False
        self.last_processing: datetime.datetime = None
        self.nomad_version: str = None
        self.nomad_commit: str = None

        # user metadata, i.e. quantities given and editable by the user
        self.with_embargo: bool = None
        self.published: bool = False
        self.coauthors: List[utils.POPO] = []
        self.shared_with: List[utils.POPO] = []
        self.comment: str = None
        self.references: List[utils.POPO] = []
        self.datasets: List[utils.POPO] = []
        self.external_id: str = None

        # parser related general (not domain specific) metadata
        self.parser_name = None

        self.update(**kwargs)

    def __getitem__(self, key):
        value = getattr(self, key, None)

        if value is None or key in ['backend']:
            raise KeyError()

        if isinstance(value, MSection):
            value = value.m_to_dict()

        return value

    def __iter__(self):
        for key, value in self.__dict__.items():
            if value is None or key in ['backend']:
                continue

            yield key

    def __len__(self):
        count = 0
        for key, value in self.__dict__.items():
            if value is None or key in ['backend']:
                continue
            count += 1

        return count

    def to_dict(self):
        return {key: value for key, value in self.items()}

    def update(self, **kwargs):
        for key, value in kwargs.items():
            if value is None:
                continue

            if isinstance(value, list):
                if len(value) == 0:
                    continue

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
        self.external_id = metadata.get('external_id')

    def apply_domain_metadata(self, backend):
        raise NotImplementedError()


class DomainQuantity:
    """
    This class can be used to define further details about a domain specific metadata
    quantity.

    Attributes:
        name: The name of the quantity, also the key used to store values in
            :class:`CalcWithMetadata`
        description: A human friendly description. The description is used to define
            the swagger documentation on the relevant API endpoints.
        multi: Indicates a list of values. This is important for the elastic mapping.
        order_default: Indicates that this metric should be used for the default order of
            search results.
        aggregations: Indicates that search aggregations (and how many) should be provided.
            0 (the default) means no aggregations.
        metric: Indicates that this quantity should be used as search metric. Values need
            to be tuples with metric name and elastic aggregation (e.g. sum, cardinality)
        zero_aggs: Return aggregation values for values with zero hits in the search. Default
            is with zero aggregations.
        elastic_mapping: An optional elasticsearch_dsl mapping. Default is ``Keyword``.
        elastic_search_type: An optional elasticsearch search type. Default is ``term``.
        elastic_field: An optional elasticsearch key. Default is the name of the quantity.
        elastic_value: A collable that takes a :class:`CalcWithMetadata` as input and produces the
            value for the elastic search index.
    """

    def __init__(
            self, description: str = None, multi: bool = False, aggregations: int = 0,
            order_default: bool = False, metric: Tuple[str, str] = None,
            zero_aggs: bool = True, metadata_field: str = None,
            elastic_mapping: type = None,
            elastic_search_type: str = 'term', elastic_field: str = None,
            elastic_value: Callable[[Any], Any] = None):

        self._name: str = None
        self.description = description
        self.multi = multi
        self.order_default = order_default
        self.aggregations = aggregations
        self.metric = metric
        self.zero_aggs = zero_aggs
        self.elastic_mapping = elastic_mapping
        self.elastic_search_type = elastic_search_type
        self.metadata_field = metadata_field
        self.elastic_field = elastic_field

        self.elastic_value = elastic_value
        if self.elastic_value is None:
            self.elastic_value = lambda o: o

        if self.elastic_mapping is None:
            self.elastic_mapping = Keyword(multi=self.multi)

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        self._name = name
        if self.metadata_field is None:
            self.metadata_field = name
        if self.elastic_field is None:
            self.elastic_field = name


class Domain:
    """
    A domain defines all metadata quantities that are specific to a certain scientific
    domain, e.g. DFT calculations, or experimental material science.

    Each domain needs to define a subclass of :class:`CalcWithMetadata`. This
    class has to define the necessary domain specific metadata quantities and how these
    are filled from parser results (usually an instance of :class:LocalBackend).

    Furthermore, the class method :func:`register_domain` of this ``Domain`` class has
    to be used to register a domain with ``domain_nam``. This also allows to provide
    further descriptions on each domain specific quantity via instance of :class:`DomainQuantity`.

    While there can be multiple domains registered. Currently, only one domain can be
    active. This active domain is define in the configuration using the ``domain_name``.

    Arguments:
        name: A name for the domain. This is used as key in the configuration ``config.domain``.
        domain_entry_class: A subclass of :class:`CalcWithMetadata` that adds the
            domain specific quantities.
        quantities: Additional specifications for the quantities in ``domain_entry_class`` as
            instances of :class:`DomainQuantity`.
        root_sections: The name of the possible root sections for this domain.
        metainfo_all_package: The name of the full metainfo package for this domain.
    """
    instance: 'Domain' = None
    instances: Dict[str, 'Domain'] = {}

    base_quantities = dict(
        authors=DomainQuantity(
            elastic_field='authors.name.keyword', multi=True, aggregations=1000,
            description=(
                'Search for the given author. Exact keyword matches in the form "Lastname, '
                'Firstname".')),
        comment=DomainQuantity(
            elastic_search_type='match', multi=True,
            description='Search within the comments. This is a text search ala google.'),
        paths=DomainQuantity(
            elastic_search_type='match', elastic_field='files', multi=True,
            description='Search for elements in one of the file paths. The paths are split at all "/".'),
        files=DomainQuantity(
            elastic_field='files.keyword', multi=True,
            description='Search for exact file name with full path.'),
        quantities=DomainQuantity(
            multi=True,
            description='Search for the existence of a certain meta-info quantity'),
        upload_id=DomainQuantity(description='Search for the upload_id.'),
        calc_id=DomainQuantity(description='Search for the calc_id.'),
        pid=DomainQuantity(description='Search for the pid.'),
        raw_id=DomainQuantity(description='Search for the raw_id.'),
        mainfile=DomainQuantity(description='Search for the mainfile.'),
        external_id=DomainQuantity(description='External user provided id. Does not have to be unique necessarily.'),
        dataset=DomainQuantity(
            elastic_field='datasets.name', multi=True, elastic_search_type='match',
            description='Search for a particular dataset by name.'),
        dataset_id=DomainQuantity(
            elastic_field='datasets.id', multi=True,
            description='Search for a particular dataset by its id.'),
        doi=DomainQuantity(
            elastic_field='datasets.doi', multi=True,
            description='Search for a particular dataset by doi (incl. http://dx.doi.org).'))

    base_metrics = dict(
        datasets=('datasets.id', 'cardinality'),
        uploaders=('uploader.name.keyword', 'cardinality'),
        authors=('authors.name.keyword', 'cardinality'),
        unique_entries=('calc_hash', 'cardinality'))

    def __init__(
            self, name: str, domain_entry_class: Type[CalcWithMetadata],
            quantities: Dict[str, DomainQuantity],
            metrics: Dict[str, Tuple[str, str]],
            default_statistics: List[str],
            root_sections=['section_run', 'section_entry_info'],
            metainfo_all_package='all.nomadmetainfo.json') -> None:

        domain_quantities = quantities
        domain_metrics = metrics

        if name == config.domain:
            assert Domain.instance is None, 'you can only define one domain.'
            Domain.instance = self

        Domain.instances[name] = self

        self.name = name
        self.domain_entry_class = domain_entry_class
        self.domain_quantities: Dict[str, DomainQuantity] = {}
        self.root_sections = root_sections
        self.metainfo_all_package = metainfo_all_package
        self.default_statistics = default_statistics

        reference_domain_calc = domain_entry_class()
        reference_general_calc = CalcWithMetadata()

        # add non specified quantities from additional metadata class fields
        for quantity_name in reference_domain_calc.__dict__.keys():
            if not hasattr(reference_general_calc, quantity_name):
                quantity = domain_quantities.get(quantity_name, None)

                if quantity is None:
                    domain_quantities[quantity_name] = DomainQuantity()

        # add all domain quantities
        for quantity_name, quantity in domain_quantities.items():
            quantity.name = quantity_name
            self.domain_quantities[quantity.name] = quantity

            # update the multi status from an example value
            if quantity.metadata_field in reference_domain_calc.__dict__:
                quantity.multi = isinstance(
                    reference_domain_calc.__dict__[quantity.metadata_field], list)

            assert not hasattr(reference_general_calc, quantity_name), \
                'quantity overrides general non domain quantity'

        # construct search quantities from base and domain quantities
        self.quantities = dict(**Domain.base_quantities)
        for quantity_name, quantity in self.quantities.items():
            quantity.name = quantity_name
        self.quantities.update(self.domain_quantities)

        assert any(quantity.order_default for quantity in Domain.instances[name].quantities.values()), \
            'you need to define a order default quantity'

        # construct metrics from base and domain metrics
        self.metrics = dict(**Domain.base_metrics)
        self.metrics.update(**domain_metrics)

    @property
    def metrics_names(self) -> Iterable[str]:
        """ Just the names of all metrics. """
        return list(self.metrics.keys())

    @property
    def aggregations(self) -> Dict[str, int]:
        """
        The search aggregations and the default maximum number of calculated buckets. See also
        :func:`nomad.search.aggregations`.
        """
        return {
            quantity.name: quantity.aggregations
            for quantity in self.quantities.values()
            if quantity.aggregations > 0
        }

    @property
    def aggregations_names(self) -> Iterable[str]:
        """ Just the names of all metrics. """
        return list(self.aggregations.keys())


def get_optional_backend_value(backend, key, section, unavailable_value=None, logger=None):
    # Section is section_system, section_symmetry, etc...
    val = None  # Initialize to None, so we can compare section values.
    # Loop over the sections with the name section in the backend.
    for section_index in backend.get_sections(section):
        if section == 'section_system':
            try:
                if not backend.get_value('is_representative', section_index):
                    continue
            except KeyError:
                continue

        try:
            new_val = backend.get_value(key, section_index)
        except KeyError:
            new_val = None

        # Compare values from iterations.
        if val is not None and new_val is not None:
            if val.__repr__() != new_val.__repr__() and logger:
                logger.warning(
                    'The values for %s differ between different %s: %s vs %s' %
                    (key, section, str(val), str(new_val)))

        val = new_val if new_val is not None else val

    if val is None and logger:
        logger.warning(
            'The values for %s where not available in any %s' % (key, section))
        return unavailable_value if unavailable_value is not None else config.services.unavailable_value
    else:
        if isinstance(val, np.generic):
            return val.item()

        return val
