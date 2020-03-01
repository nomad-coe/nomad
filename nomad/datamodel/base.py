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

import numpy as np

from nomad import config

# from .metainfo import Dataset, User, EntryMetadata


# class DomainQuantity:
#     '''
#     This class can be used to define further details about a domain specific metadata
#     quantity.

#     Attributes:
#         name: The name of the quantity, also the key used to store values in
#             :class:`EntryMetadata`
#         description: A human friendly description. The description is used to define
#             the swagger documentation on the relevant API endpoints.
#         multi: Indicates a list of values. This is important for the elastic mapping.
#         order_default: Indicates that this metric should be used for the default order of
#             search results.
#         aggregations: Indicates that search aggregations (and how many) should be provided.
#             0 (the default) means no aggregations.
#         metric: Indicates that this quantity should be used as search metric. Values need
#             to be tuples with metric name and elastic aggregation (e.g. sum, cardinality)
#         elastic_mapping: An optional elasticsearch_dsl mapping. Default is ``Keyword``.
#         elastic_search_type: An optional elasticsearch search type. Default is ``term``.
#         elastic_field: An optional elasticsearch key. Default is the name of the quantity.
#         elastic_value: A collable that takes a :class:`EntryMetadata` as input and produces the
#             value for the elastic search index.
#         argparse_action: Action to use on argparse, either append or split for multi values. Append is default.
#     '''

#     def __init__(
#             self, description: str = None, multi: bool = False, aggregations: int = 0,
#             order_default: bool = False, metric: Tuple[str, str] = None,
#             metadata_field: str = None, elastic_mapping: type = None,
#             elastic_search_type: str = 'term', elastic_field: str = None,
#             elastic_value: Callable[[Any], Any] = None,
#             argparse_action: str = 'append'):

#         self.domain: str = None
#         self._name: str = None
#         self.description = description
#         self.multi = multi
#         self.order_default = order_default
#         self.aggregations = aggregations
#         self.metric = metric
#         self.elastic_mapping = elastic_mapping
#         self.elastic_search_type = elastic_search_type
#         self.metadata_field = metadata_field
#         self.elastic_field = elastic_field
#         self.argparse_action = argparse_action

#         self.elastic_value = elastic_value
#         if self.elastic_value is None:
#             self.elastic_value = lambda o: o

#         if self.elastic_mapping is None:
#             self.elastic_mapping = Keyword(multi=self.multi)

#     @property
#     def name(self) -> str:
#         return self._name

#     @name.setter
#     def name(self, name: str) -> None:
#         self._name = name
#         if self.metadata_field is None:
#             self.metadata_field = name
#         if self.elastic_field is None:
#             self.elastic_field = self.name

#     @property
#     def qualified_elastic_field(self) -> str:
#         if self.domain is None:
#             return self.elastic_field
#         else:
#             return '%s.%s' % (self.domain, self.elastic_field)

#     @property
#     def qualified_name(self) -> str:
#         if self.domain is None:
#             return self.name
#         else:
#             return '%s.%s' % (self.domain, self.name)


# def only_atoms(atoms):
#     numbers = [ase.data.atomic_numbers[atom] for atom in atoms]
#     only_atoms = [ase.data.chemical_symbols[number] for number in sorted(numbers)]
#     return ''.join(only_atoms)


# class Domain:
#     '''
#     A domain defines all metadata quantities that are specific to a certain scientific
#     domain, e.g. DFT calculations, or experimental material science.

#     Each domain needs to define a subclass of :class:`EntryMetadata`. This
#     class has to define the necessary domain specific metadata quantities and how these
#     are filled from parser results (usually an instance of :class:LocalBackend).

#     Furthermore, the class method :func:`register_domain` of this ``Domain`` class has
#     to be used to register a domain with ``domain_nam``. This also allows to provide
#     further descriptions on each domain specific quantity via instance of :class:`DomainQuantity`.

#     While there can be multiple domains registered. Currently, only one domain can be
#     active. This active domain is define in the configuration using the ``domain_name``.

#     Arguments:
#         name: A name for the domain. This is used as key in the configuration ``config.domain``.
#         domain_entry_class: A subclass of :class:`EntryMetadata` that adds the
#             domain specific quantities.
#         quantities: Additional specifications for the quantities in ``domain_entry_class`` as
#             instances of :class:`DomainQuantity`.
#         metrics: Tuples of elastic field name and elastic aggregation operation that
#             can be used to create statistic values.
#         groups: Tuple of quantity name and metric that describes quantities that
#             can be used to group entries by quantity values.
#         root_sections: The name of the possible root sections for this domain.
#         metainfo_all_package: The name of the full metainfo package for this domain.
#     '''
#     instances: Dict[str, 'Domain'] = {}

#     base_quantities = dict(
#         authors=DomainQuantity(
#             elastic_field='authors.name.keyword', multi=True, aggregations=1000,
#             description=(
#                 'Search for the given author. Exact keyword matches in the form "Lastname, '
#                 'Firstname".')),
#         uploader_id=DomainQuantity(
#             elastic_field='uploader.user_id', multi=False, aggregations=5,
#             description=('Search for the given uploader id.')),
#         uploader_name=DomainQuantity(
#             elastic_field='uploader.name.keyword', multi=False,
#             description=('Search for the exact uploader\'s full name')),
#         comment=DomainQuantity(
#             elastic_search_type='match', multi=True,
#             description='Search within the comments. This is a text search ala google.'),
#         paths=DomainQuantity(
#             elastic_search_type='match', elastic_field='files', multi=True,
#             description='Search for elements in one of the file paths. The paths are split at all "/".'),
#         files=DomainQuantity(
#             elastic_field='files.keyword', multi=True,
#             description='Search for exact file name with full path.'),
#         quantities=DomainQuantity(
#             multi=True,
#             description='Search for the existence of a certain meta-info quantity'),
#         upload_id=DomainQuantity(
#             description='Search for the upload_id.',
#             multi=True, argparse_action='split', elastic_search_type='terms'),
#         upload_time=DomainQuantity(
#             description='Search for the exact upload time.', elastic_search_type='terms'),
#         upload_name=DomainQuantity(
#             description='Search for the upload_name.',
#             multi=True, argparse_action='split', elastic_search_type='terms'),
#         calc_id=DomainQuantity(
#             description='Search for the calc_id.',
#             multi=True, argparse_action='split', elastic_search_type='terms'),
#         pid=DomainQuantity(
#             description='Search for the pid.',
#             multi=True, argparse_action='split', elastic_search_type='terms'),
#         raw_id=DomainQuantity(
#             description='Search for the raw_id.',
#             multi=True, argparse_action='split', elastic_search_type='terms'),
#         mainfile=DomainQuantity(
#             description='Search for the mainfile.',
#             multi=True, argparse_action='append', elastic_search_type='terms'),
#         external_id=DomainQuantity(
#             description='External user provided id. Does not have to be unique necessarily.',
#             multi=True, argparse_action='split', elastic_search_type='terms'),
#         calc_hash=DomainQuantity(
#             description='Search for the entries hash.',
#             multi=True, argparse_action='split', elastic_search_type='terms'),
#         dataset=DomainQuantity(
#             elastic_field='datasets.name', multi=True, elastic_search_type='match',
#             description='Search for a particular dataset by name.'),
#         dataset_id=DomainQuantity(
#             elastic_field='datasets.id', multi=True,
#             description='Search for a particular dataset by its id.'),
#         doi=DomainQuantity(
#             elastic_field='datasets.doi', multi=True,
#             description='Search for a particular dataset by doi (incl. http://dx.doi.org).'),
#         formula=DomainQuantity(
#             'The chemical (hill) formula of the simulated system.',
#             order_default=True),
#         atoms=DomainQuantity(
#             'The atom labels of all atoms in the simulated system.',
#             aggregations=len(ase.data.chemical_symbols), multi=True),
#         only_atoms=DomainQuantity(
#             'The atom labels concatenated in species-number order. Used with keyword search '
#             'to facilitate exclusive searches.',
#             elastic_value=only_atoms, metadata_field='atoms', multi=True),
#         n_atoms=DomainQuantity(
#             'Number of atoms in the simulated system',
#             elastic_mapping=Integer()))

#     base_metrics = dict(
#         datasets=('dataset_id', 'cardinality'),
#         uploads=('upload_id', 'cardinality'),
#         uploaders=('uploader_name', 'cardinality'),
#         authors=('authors', 'cardinality'),
#         unique_entries=('calc_hash', 'cardinality'))

#     base_groups = dict(
#         datasets=('dataset_id', 'datasets'),
#         uploads=('upload_id', 'uploads'))

#     @classmethod
#     def get_quantity(cls, name_spec) -> DomainQuantity:
#         '''
#         Returns the quantity definition for the given quantity name. The name can be the
#         qualified name (``domain.quantity``) or in Django-style (``domain__quantity``).
#         '''
#         qualified_name = name_spec.replace('__', '.')
#         split_name = qualified_name.split('.')
#         if len(split_name) == 1:
#             return cls.base_quantities[split_name[0]]
#         elif len(split_name) == 2:
#             return cls.instances[split_name[0]].quantities[split_name[1]]
#         else:
#             assert False, 'qualified quantity name depth must be 2 max'

#     @classmethod
#     def all_quantities(cls) -> Iterable[DomainQuantity]:
#         return set([quantity for domain in cls.instances.values() for quantity in domain.quantities.values()])

#     def __init__(
#             self, name: str, domain_entry_class: Type[EntryMetadata],
#             quantities: Dict[str, DomainQuantity],
#             metrics: Dict[str, Tuple[str, str]],
#             groups: Dict[str, Tuple[str, str]],
#             default_statistics: List[str],
#             root_sections=['section_run', 'section_entry_info'],
#             metainfo_all_package='all.nomadmetainfo.json') -> None:

#         domain_quantities = quantities

#         Domain.instances[name] = self

#         self.name = name
#         self.domain_entry_class = domain_entry_class
#         self.domain_quantities: Dict[str, DomainQuantity] = {}
#         self.root_sections = root_sections
#         self.metainfo_all_package = metainfo_all_package
#         self.default_statistics = default_statistics

#         # TODO
#         return

#         reference_domain_calc = EntryMetadata(domain=name)
#         reference_general_calc = EntryMetadata(domain=None)

#         # add non specified quantities from additional metadata class fields
#         for quantity_name in reference_domain_calc.__dict__.keys():
#             if not hasattr(reference_general_calc, quantity_name):
#                 quantity = domain_quantities.get(quantity_name, None)
#                 if quantity is None:
#                     domain_quantities[quantity_name] = DomainQuantity()

#         # ensure domain quantity names and domains
#         for quantity_name, quantity in domain_quantities.items():
#             quantity.domain = name
#             quantity.name = quantity_name

#         # add domain prefix to domain metrics and groups
#         domain_metrics = {
#             '%s.%s' % (name, key): (quantities[quantity].qualified_elastic_field, es_op)
#             for key, (quantity, es_op) in metrics.items()}
#         domain_groups = {
#             '%s.%s' % (name, key): (quantities[quantity].qualified_name, '%s.%s' % (name, metric))
#             for key, (quantity, metric) in groups.items()}

#         # add all domain quantities
#         for quantity_name, quantity in domain_quantities.items():
#             self.domain_quantities[quantity.name] = quantity

#             # update the multi status from an example value
#             if quantity.metadata_field in reference_domain_calc.__dict__:
#                 quantity.multi = isinstance(
#                     reference_domain_calc.__dict__[quantity.metadata_field], list)

#             assert not hasattr(reference_general_calc, quantity_name), \
#                 'quantity overrides general non domain quantity: %s' % quantity_name

#         # construct search quantities from base and domain quantities
#         self.quantities = dict(**Domain.base_quantities)
#         for quantity_name, quantity in self.quantities.items():
#             quantity.name = quantity_name
#         self.quantities.update(self.domain_quantities)

#         assert any(quantity.order_default for quantity in Domain.instances[name].quantities.values()), \
#             'you need to define a order default quantity'

#         # construct metrics from base and domain metrics
#         self.metrics = dict(**Domain.base_metrics)
#         self.metrics.update(**domain_metrics)
#         self.groups = dict(**Domain.base_groups)
#         self.groups.update(**domain_groups)

#     @property
#     def metrics_names(self) -> Iterable[str]:
#         ''' Just the names of all metrics. '''
#         return list(self.metrics.keys())

#     @property
#     def aggregations(self) -> Dict[str, int]:
#         '''
#         The search aggregations and the default maximum number of calculated buckets. See also
#         :func:`nomad.search.aggregations`.
#         '''
#         return {
#             quantity.name: quantity.aggregations
#             for quantity in self.quantities.values()
#             if quantity.aggregations > 0
#         }

#     @property
#     def aggregations_names(self) -> Iterable[str]:
#         ''' Just the names of all metrics. '''
#         return list(self.aggregations.keys())

#     @property
#     def order_default_quantity(self) -> str:
#         for quantity in self.quantities.values():
#             if quantity.order_default:
#                 return quantity.qualified_name

#         assert False, 'each domain must defina an order_default quantity'


def get_optional_backend_value(backend, key, section, unavailable_value=None, logger=None):
    # Section is section_system, section_symmetry, etc...
    val = None  # Initialize to None, so we can compare section values.
    # Loop over the sections with the name section in the backend.
    for section_index in backend.get_sections(section):
        if section == 'section_system':
            try:
                if not backend.get_value('is_representative', section_index):
                    continue
            except (KeyError, IndexError):
                continue

        try:
            new_val = backend.get_value(key, section_index)
        except (KeyError, IndexError):
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
