from nomad.bundles import get_section_defs_for_upload
from nomad.datamodel import EntryArchive
from nomad.datamodel.datamodel import EntryMetadata


def test_get_section_defs_for_upload(non_empty_processed_with_temporal):
    indexed_definition_ids = set()
    entry_ids = [
        entry.entry_id for entry in non_empty_processed_with_temporal.successful_entries
    ]

    with non_empty_processed_with_temporal.entries_metadata(
        entry_ids=entry_ids
    ) as entries_metadata:
        for entry_metadata in entries_metadata:
            for section_def in entry_metadata.section_defs:
                indexed_definition_ids.add(section_def.definition_id)

    definitions = get_section_defs_for_upload(
        non_empty_processed_with_temporal.upload_id
    )

    definition_ids = [definition.definition_id for definition in definitions]
    assert len(definition_ids) == len(set(definition_ids))
    assert set(definition_ids) == indexed_definition_ids

    qualified_names = {definition.qualified_name() for definition in definitions}
    assert EntryArchive.m_def.qualified_name() in qualified_names
    assert EntryMetadata.m_def.qualified_name() in qualified_names
    assert 'nomad.datamodel.results.Results' in qualified_names
    assert 'runschema.run.Run' in qualified_names
