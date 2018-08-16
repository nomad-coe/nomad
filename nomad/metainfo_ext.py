repository_data = [
    {
        'name': 'Code',
        'metainfo_name': 'program_name',
        'metainfo_path': ['section_run'],
        'search': True,
        'data_type': str,
    },
    {
        'name': 'Version',
        'metainfo_name': 'program_version',
        'metainfo_path': ['section_run'],
        'comments': 'Should be treated via regexp.',
        'data_type': str,
    },
    {

    }

    # method
    # program_basis_set_type
    # XC_functional_name*
    # electronic_structure_method

]

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
# IndexEntry("system_type"),            // from normaliser system-type  +++
# IndexEntry("crystal_system"),         // from normaliser symmetry     +++
# IndexEntry("space_group_number"),     // from normaliser symmetry     +++
# IndexEntry("springer_id"),            // from normaliser springer
# IndexEntry("springer_classification"),// from normaliser springer
# IndexEntry("configuration_raw_gid"),  // from normaliser stats            DONE
# //computational setup
# IndexEntry("electronic_structure_method"),                            +++
# IndexEntry("XC_functional_name"),     // from normaliser stats does XC_functional_type  +++
# IndexEntry("program_basis_set_type"),                                 +++
# //data related to other projects ingested into Nomad(e.g. AFlowLib)
# IndexEntry("prototype_aflow_id"),     // from normalizer prototypes
# IndexEntry("prototype_aflow_url"),    // from normalizer prototypes
# //auxiliary data of the query
# IndexEntry("number_of_single_configurations", Some(intField), Some { calc =>
#     Seq(calc.sectionTable(Seq("section_run", "section_system")).lengthL)
# })


# val activeNormalizers: Map[String, NormalizerGenerator] = {
#     listToMap(Seq(
#       StatsNormalizer,            # python
#       FhiAimsBasisNormalizer,     # python
#       SymmetryNormalizer,         # python
#       SpringerNormalizer,         # python
#       MetaIndexRebuildNormalizer, # scala
#       SectionSystemNormalizer,    # scala
#       DosValuesByAtomNormalizer,  # scala
#       RepoTagsNormalizer,         # python
#       PrototypesNormalizer,       # python
#       BandStructureNormalizer,    # python
#       FhiDosNormalizer,           # scala
#       CalculationInfoNormalizer,  # scala
#       SystemTypeNormalizer,       # python
#       UploadInfoNormalizer        # scala
#     ))
#   }