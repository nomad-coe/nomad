from elasticsearch_dsl import Document, Date, Keyword, connections
import logging
import inspect

from nomad.parsing import LocalBackend

logger = logging.getLogger(__name__)
connections.create_connection(hosts=['localhost'])


class Calc(Document):
    """
    The search index representation of a calculation. All information for this should
    be available after parsing and normalization. It should contain all fields and
    mappings that are necessary to create the desired repository/archive functionality.
    """
    calc_hash = Keyword()

    upload_time = Date()
    upload_hash = Keyword()
    mainfile = Keyword()

    program_name = Keyword()
    program_version = Keyword()

    atom_species = Keyword()
    system_type = Keyword()
    crystal_system = Keyword()
    space_group_number = Keyword()
    configuration_raw_gid = Keyword()
    XC_functional_name = Keyword()

    class Index:
        name = 'calcs'

    @staticmethod
    def add_from_backend(backend: LocalBackend, **kwargs) -> 'Calc':
        """
        Add the calc data from the given backend to the elastic search index. Additional
        meta-data can be given as *kwargs*.
        """

        upload_hash = kwargs.get('upload_hash', None)
        calc_hash = kwargs.get('calc_hash', None)

        assert upload_hash is not None and calc_hash is not None

        calc = Calc(meta=dict(id='%s/%s' % (upload_hash, calc_hash)))

        for property in Calc._doc_type.mapping:
            if property in kwargs:
                value = kwargs[property]
            else:
                try:
                    value = backend.get_value(property, 0)
                except KeyError:
                    logger.warning(
                        'No value for property %s could be extracted for calc %s/%s.' %
                        (property, upload_hash, calc_hash))
                    continue

            setattr(calc, property, value)

        calc.save()
        return calc


Calc.init()

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