from typing import Generator
import pytest
from minio import ResponseError

import nomad.files as files
import nomad.config as config
from nomad.processing import UploadProcessing

example_upload_id = '__test_upload_id'
example_file = 'data/examples_vasp.zip'


@pytest.fixture
def uploaded_id() -> Generator[str, None, None]:
    files._client.fput_object(config.s3.uploads_bucket, example_upload_id, example_file)

    yield example_upload_id

    try:
        files._client.remove_object(config.s3.uploads_bucket, example_upload_id)
    except ResponseError:
        pass


def test_processing(uploaded_id):
        run = UploadProcessing(uploaded_id)
        run.start()
        run.get(timeout=30)
        assert run.ready()

        run.forget()

        assert run.task_name == 'nomad.processing.close_upload'
        assert run.status == 'SUCCESS'
