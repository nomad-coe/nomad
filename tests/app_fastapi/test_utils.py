from typing import Iterator
import os.path
import zipfile

from nomad import config
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.app_fastapi.utils import create_streamed_zipfile, File

from tests.conftest import clear_raw_files
from tests.test_files import create_test_upload_files


def test_create_streamed_zip(raw_files_infra):
    # We use the files of a simpe test upload to create streamed zip with all the raw
    # files.
    archive = EntryArchive()
    metadata = archive.m_create(EntryMetadata)
    metadata.upload_id = 'test_id'
    metadata.calc_id = 'test_id'
    metadata.mainfile = 'root/subdir/mainfile.json'

    upload_files = create_test_upload_files('test_id', [archive])

    def generate_files() -> Iterator[File]:
        for path in upload_files.raw_file_manifest():
            with upload_files.raw_file(path) as f:
                yield File(
                    path=path,
                    f=f,
                    size=upload_files.raw_file_size(path))

    if not os.path.exists(config.fs.tmp):
        os.makedirs(config.fs.tmp)

    zip_file_path = os.path.join(config.fs.tmp, 'results.zip')
    with open(zip_file_path, 'wb') as f:
        for content in create_streamed_zipfile(generate_files()):
            f.write(content)

    with zipfile.ZipFile(zip_file_path) as zf:
        assert zf.testzip() is None

    clear_raw_files()
