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
The interface towards our search engine. It allows to store calculations as documents
of search relevant properties.

..autoclass:: nomad.search.Calc:
"""

import sys
from datetime import datetime
import elasticsearch.exceptions
from elasticsearch_dsl import Document, Date, Keyword, Search, connections

from nomad import config, files
from nomad.parsing import LocalBackend
from nomad.utils import get_logger

logger = get_logger(__name__)

# ensure elastic connection
if 'sphinx' not in sys.modules:
    client = connections.create_connection(hosts=[config.elastic.host])


key_mappings = {
    'basis_set_type': 'program_basis_set_type',
    'chemical_composition': 'chemical_composition_bulk_reduced'
}


class AlreadyExists(Exception): pass


class Calc(Document):
    """
    Instances of this class represent calculations. This class manages the elastic
    search index entry, files, and archive for the respective calculation.
    Instances should be created directly, but via the static :func:`create_from_backend`
    function.

    The attribute list, does not include the various repository properties generated
    while parsing, including ``program_name``, ``program_version``, etc.

    Attributes:
        calc_hash: The hash that identified the calculation within an upload
        upload_hash: The hash of the upload
        upload_id: The id of the upload used to create this calculation
        mainfile: The mainfile (including path in upload) that was used to create this calc
    """
    calc_hash = Keyword()

    upload_time = Date()
    upload_id = Keyword()
    upload_hash = Keyword()
    mainfile = Keyword()

    program_name = Keyword()
    program_version = Keyword()

    chemical_composition = Keyword()
    basis_set_type = Keyword()
    atom_species = Keyword()
    system_type = Keyword()
    crystal_system = Keyword()
    space_group_number = Keyword()
    configuration_raw_gid = Keyword()
    XC_functional_name = Keyword()

    class Index:
        name = config.elastic.calc_index

    @property
    def archive_id(self) -> str:
        """ The unique id for this calculation. """
        return '%s/%s' % (self.upload_hash, self.calc_hash)

    def delete(self):
        """
        Delete this calculation and all associated data. This includes all files,
        the archive, and this search index entry.
        """
        # delete the archive
        files.delete_archive(self.archive_id)

        # delete the search index entry
        super().delete()

    @staticmethod
    def es_search(body):
        """ Perform an elasticsearch and not elasticsearch_dsl search on the Calc index. """
        return client.search(index=config.elastic.calc_index, body=body)

    @staticmethod
    def delete_all(**kwargs):
        for calc in Calc.search().query('match', **kwargs).execute():
            calc.delete()

    @staticmethod
    def create_from_backend(
            backend: LocalBackend, upload_id: str, upload_hash: str, calc_hash: str, **kwargs) \
            -> 'Calc':
        """
        Create a new calculation instance. The data from the given backend
        will be used. Additional meta-data can be given as *kwargs*. ``upload_id``,
        ``upload_hash``, and ``calc_hash`` are mandatory.
        This will create a elastic search entry and store the backend data to the
        archive.

        Arguments:
            backend: The parsing/normalizing backend that contains the calculation data.
            upload_hash: The upload hash of the originating upload.
            upload_id: The upload id of the originating upload.
            calc_hash: The upload unique hash for this calculation.
            kwargs: Additional arguments not stored in the backend.

        Raises:
            AlreadyExists: If the calculation already exists in elastic search. We use
                the elastic document lock here. The elastic document is ided via the
                ``archive_id``.
        """
        assert upload_hash is not None and calc_hash is not None and upload_id is not None
        kwargs.update(dict(upload_hash=upload_hash, calc_hash=calc_hash, upload_id=upload_id))

        calc = Calc(meta=dict(id='%s/%s' % (upload_hash, calc_hash)))

        for property in Calc._doc_type.mapping:
            property = key_mappings.get(property, property)

            if property in kwargs:
                value = kwargs[property]
            else:
                try:
                    value = backend.get_value(property, 0)
                except KeyError:
                    logger.warning(
                        'Missing property value', property=property, upload_id=upload_id,
                        upload_hash=upload_hash, calc_hash=calc_hash)
                    continue

            setattr(calc, property, value)

        # persist to elastic search
        try:
            calc.save(op_type='create')
        except Exception as e:
            raise AlreadyExists('Calculation %s does already exist.' % (calc.archive_id))

        # persist the archive
        with files.write_archive_json(calc.archive_id) as out:
            backend.write_json(out, pretty=True)

        return calc

    @property
    def json_dict(self):
        """ A json serializable dictionary representation. """
        data = self.to_dict()

        upload_time = data.get('upload_time', None)
        if upload_time is not None and isinstance(upload_time, datetime):
            data['upload_time'] = data['upload_time'].isoformat()

        data['archive_id'] = self.archive_id

        return {key: value for key, value in data.items() if value is not None}

    @staticmethod
    def upload_exists(upload_hash):
        """ Returns true if there are already calcs from the given upload. """
        search = Search(using=client, index=config.elastic.calc_index) \
            .query('match', upload_hash=upload_hash) \
            .execute()

        return len(search) > 0


if 'sphinx' not in sys.modules:
    try:
        Calc.init()
    except elasticsearch.exceptions.RequestError as e:
        if e.status_code == 400 and 'resource_already_exists_exception' in e.error:
            pass  # happens if two services try this at the same time
        else:
            raise e


# Taken from the IndexManifest (scala, NOMAD-coe)
#
# //calculation identifiers (these are treated separately in the importer)
# IndexEntry("calculation_gid"),                                        +++
# IndexEntry("archive_gid"),                                            +++
# //calculation data
# IndexEntry("main_file_uri"),                                          +++
# IndexEntry("calculation_uploader_name"),
# IndexEntry("calculation_upload_date"),
# IndexEntry("calculation_pid"),
# IndexEntry("program_name"),                                           +++ + version
# IndexEntry("stats_meta_present"),
# //system related data
# IndexEntry("atom_species"),           // from normaliser stats        +++ DONE
# IndexEntry("system_composition"),     // from normaliser springer
# IndexEntry("system_reweighted_composition"),  // ???
# IndexEntry("system_type"),            // from normaliser system-type  +++ DONE
# IndexEntry("crystal_system"),         // from normaliser symmetry     +++ DONE
# IndexEntry("space_group_number"),     // from normaliser symmetry     +++ DONE
# IndexEntry("springer_id"),            // from normaliser springer
# IndexEntry("springer_classification"),// from normaliser springer
# IndexEntry("configuration_raw_gid"),  // from normaliser stats            DONE
# //computational setup
# IndexEntry("electronic_structure_method"),                            +++
# IndexEntry("XC_functional_name"),     // from normaliser stats does XC_functional_type  +++ DONE
# IndexEntry("program_basis_set_type"),                                 +++
# //data related to other projects ingested into Nomad(e.g. AFlowLib)
# IndexEntry("prototype_aflow_id"),     // from normalizer prototypes
# IndexEntry("prototype_aflow_url"),    // from normalizer prototypes
# //auxiliary data of the query
# IndexEntry("number_of_single_configurations", Some(intField), Some { calc =>
#     Seq(calc.sectionTable(Seq("section_run", "section_system")).lengthL)
# })
