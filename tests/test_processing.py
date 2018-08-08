import unittest
from unittest import TestCase
import time
import logging
from minio import ResponseError

import nomad.files as files
import nomad.config as config
from nomad.processing import start_process_upload, get_process_upload_state

test_upload_id = '__test_upload_id'


class ProcessingTests(TestCase):

    def setUp(self):
        files._client.fput_object(config.s3.uploads_bucket, test_upload_id, 'data/examples_vasp.zip')

    def tearDown(self):
        try:
            files._client.remove_object(config.s3.uploads_bucket, test_upload_id)
        except ResponseError:
            pass

    def test_processing(self):
        task = start_process_upload(test_upload_id)

        result = None
        while(True):
            time.sleep(0.0001)
            new_result = get_process_upload_state(task)
            if result != new_result:
                result = new_result
                if result['close'] == 'SUCCESS' or result['close'] == 'FAILURE':
                    break

        self.assertTrue(result['close'] == 'SUCCESS')


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()
