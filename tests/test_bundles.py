#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import io
import os
import zipfile

from nomad.bundles import BundleExporter, BundleImporter, get_section_defs_for_upload
from nomad.config import config
from nomad.datamodel import EntryArchive
from nomad.datamodel.datamodel import EntryMetadata
from nomad.processing import Entry, Upload


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


# Test bundle export


def test_archive_preserving_bundle_export_as_stream(non_empty_processed_with_temporal):
    """Test for `include_archive_files=True`"""
    # Create (export) the bundle
    exporter = BundleExporter(
        upload=non_empty_processed_with_temporal,
        export_as_stream=True,
        export_path=None,
        zipped=True,
        overwrite=False,  # not applicable for streaming
        export_settings=config.bundle_export.default_settings,
    )
    assert exporter.export_settings.include_archive_files is True

    bundle_stream = exporter.export_bundle()
    assert bundle_stream is not None

    with zipfile.ZipFile(io.BytesIO(b''.join(bundle_stream))) as zf:
        names = set(zf.namelist())

        expected_stable = {
            'bundle_info.json',
            'raw/examples_template/0.aux',
            'raw/examples_template/1.aux',
            'raw/examples_template/2.aux',
            'raw/examples_template/3.aux',
            'raw/examples_template/template.json',
        }
        assert expected_stable.issubset(names)

        archive_files = {name for name in names if name.startswith('archive/')}
        assert len(archive_files) == 1


# Tests for bundle transfer roundtrip


def test_archive_preserving_bundle_roundtrip(
    non_empty_processed_with_temporal, tmp_path
):
    """A bundle with archive files should import without reprocessing."""

    # Export bundle to file
    exporter = BundleExporter(
        upload=non_empty_processed_with_temporal,
        export_as_stream=True,
        export_path=None,
        zipped=True,
        overwrite=False,
        export_settings=config.bundle_export.default_settings,
    )
    bundle_bytes = b''.join(exporter.export_bundle())

    bundle_zip_path = tmp_path / 'bundle.zip'
    bundle_zip_path.write_bytes(bundle_bytes)

    extracted_bundle_path = tmp_path / 'bundle'
    with zipfile.ZipFile(bundle_zip_path) as zf:
        zf.extractall(extracted_bundle_path)

    non_empty_processed_with_temporal.delete_upload_local()

    # Import bundle and check
    importer = BundleImporter(
        None,
        config.bundle_import.default_settings.customize(
            dict(
                trigger_processing=False,
                delete_bundle_on_success=False,
                delete_bundle_on_fail=False,
            )
        ),
    )
    importer.open(str(extracted_bundle_path))
    try:
        imported_upload = importer.create_upload_skeleton()
        importer.import_bundle(imported_upload, True)
    finally:
        importer.close()

    imported_upload = Upload.get(imported_upload.upload_id)
    assert imported_upload.upload_id == non_empty_processed_with_temporal.upload_id

    imported_entries = list(Entry.objects(upload_id=imported_upload.upload_id))
    assert len(imported_entries) == len(
        non_empty_processed_with_temporal.successful_entries
    )

    imported_files = set()
    for dirpath, _, filenames in os.walk(imported_upload.upload_files.os_path):
        for filename in filenames:
            imported_files.add(
                os.path.relpath(
                    os.path.join(dirpath, filename),
                    imported_upload.upload_files.os_path,
                )
            )

    expected_stable = {
        'raw/examples_template/0.aux',
        'raw/examples_template/1.aux',
        'raw/examples_template/2.aux',
        'raw/examples_template/3.aux',
        'raw/examples_template/template.json',
    }
    assert expected_stable.issubset(imported_files)

    archive_files = {name for name in imported_files if name.startswith('archive/')}
    assert len(archive_files) > 0
