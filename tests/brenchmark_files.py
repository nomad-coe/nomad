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

import time

from nomad import files


def read_all(upload_id: str = 'oHbU2EjwR8Gqw8H1EE3xeQ', prefix: str = 'aflowlib_data/LIB1_LIB/Ac_s:PAW_GGA:11Apr2000/A2/'):
    upload_files = files.UploadFiles.get(upload_id)

    n_bytes = 0
    n_files = 0
    start = time.time()
    for path in upload_files.raw_file_manifest(prefix):
        try:
            with upload_files.raw_file(path, 'rb') as f:
                n_bytes += len(f.read())
                n_files += 1
        except files.Restricted:
            pass

        print('%d files, %d bytes, %f bytes/s' % (n_files, n_bytes, n_bytes / (time.time() - start)))
