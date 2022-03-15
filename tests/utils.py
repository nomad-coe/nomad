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

''' Methods to help with testing of nomad@FAIRDI.'''

from typing import List, Union, Dict, Any
import urllib.parse
import json
from logging import LogRecord
import zipfile
import os.path


def assert_log(caplog, level: str, event_part: str) -> LogRecord:
    '''
    Assert whether a log message exists in the logs of the tests at a certain level.

    Parameters
    ----------
    caplog : pytest fixture
        This informs pytest that we want to access the logs from a pytest test.
    level : str
        The level of type of log for which we will search (e.g. 'WARN',
        'ERROR', 'DEBUG').
    event_part : str
        The error message we're after. We search the logs matching level if they
        contain this string.

    '''
    record = None
    for record in caplog.get_records(when='call'):
        if record.levelname == level:
            try:
                event_data = json.loads(record.msg)
                present = event_part in event_data['event']
            except Exception:
                present = event_part in record.msg

            if present:
                record = record
                # No need to look for more matches since we aren't counting matches.
                break
    assert record is not None

    return record


def assert_at_least(source, target):
    '''
    Compares two dicts recursively and asserts that all information in source equals
    the same information in target. Additional information in target is ignored.
    '''
    for key, value in source.items():
        assert key in target, '%s with value %s in %s is not in %s' % (key, source[key], source, target)
        if isinstance(value, dict):
            assert_at_least(value, target[key])
        else:
            assert value == target[key], '%s with value %s in %s is not equal the target value %s in %s' % (
                key, source[key], source, target[key], target)


def assert_url_query_args(url: str, **kwargs):
    '''
    Parses the url, and checks that the query arguments match the values specified by kwargs.
    '''
    __, __, __, __, query, __ = urllib.parse.urlparse(url)
    query_dict = urllib.parse.parse_qs(query)
    for k, v in kwargs.items():
        if v is None:
            assert k not in query_dict
        else:
            assert query_dict[k][0] == str(v)


def build_url(base_url: str, query_args: Dict[str, Any]) -> str:
    '''
    Takes a base_url and a dictionary, and combines to a url with query arguments.
    Arguments with value None are ignored.
    '''
    # Remove args with value None
    query_args_clean = {k: v for k, v in query_args.items() if v is not None}
    if not query_args_clean:
        return base_url
    return base_url + '?' + urllib.parse.urlencode(query_args_clean, doseq=True)


def set_upload_entry_metadata(upload, metadata: Dict[str, Any]):
    '''
    Sets the provided metadata values on all entries of the given upload.
    '''
    from nomad import processing as proc
    for entry in proc.Entry.objects(upload_id=upload.upload_id):
        entry.set_mongo_entry_metadata(**metadata)
        entry.save()


def create_template_upload_file(
        tmp, mainfiles: Union[str, List[str]] = None, auxfiles: int = 4,
        directory: str = 'examples_template', name: str = 'examples_template.zip',
        more_files: Union[str, List[str]] = None):

    '''
    Creates a temporary upload.zip file based on template.json (for the artificial test
    parser) that can be used for test processings.
    '''

    if mainfiles is None:
        mainfiles = 'tests/data/proc/templates/template.json'

    if isinstance(mainfiles, str):
        mainfiles = [mainfiles]

    if more_files is None:
        more_files = []

    if isinstance(more_files, str):
        more_files = [more_files]

    upload_path = os.path.join(tmp, name)
    with zipfile.ZipFile(upload_path, 'w') as zf:
        for i in range(0, auxfiles):
            with zf.open(f'{directory}/{i}.aux', 'w') as f:
                f.write(b'content')
            for mainfile in mainfiles:
                zf.write(mainfile, f'{directory}/{os.path.basename(mainfile)}')

        for additional_file in more_files:
            zf.write(additional_file, f'{directory}/{os.path.basename(additional_file)}')

    return upload_path
