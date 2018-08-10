import unittest
from unittest import TestCase
import logging
from minio import ResponseError
import os

import nomad.files as files
import nomad.config as config
from nomad.processing import ProcessRun

test_upload_id = '__test_upload_id'


class ProcessingTests(TestCase):

    def setUp(self):
        example_file = 'data/examples_vasp.zip'
        with open(example_file, 'rb') as file_data:
            file_stat = os.stat(example_file)
            files._client.put_object(
                config.s3.uploads_bucket, test_upload_id, file_data, file_stat.st_size)

    def tearDown(self):
        try:
            files._client.remove_object(config.s3.uploads_bucket, test_upload_id)
        except ResponseError:
            pass

    def test_processing(self):
        run = ProcessRun(test_upload_id)
        run.start()
        state = run.get(timeout=30)
        self.assertTrue(run.ready())
        run.forget()

        self.assertEqual('nomad.processing.close_upload', state.task_name)
        self.assertEqual('SUCCESS', state.status)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    unittest.main()
