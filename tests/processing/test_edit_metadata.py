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
import pytest
from datetime import datetime

from nomad import datamodel, metainfo
from nomad.processing import Upload
from nomad.processing.data import _editable_metadata
from nomad.search import search
from nomad.app.v1.models import MetadataEditRequest


def assert_edit_request(user, **kwargs):
    # Extract test parameters (lots of defaults)
    upload_id = kwargs.get('upload_id', 'id_unpublished_w')
    query = kwargs.get('query')
    metadata = kwargs.get('metadata')
    entries = kwargs.get('entries')
    verify = kwargs.get('verify', False)
    expected_status_code = kwargs.get('expected_status_code', 200)
    if expected_status_code == 200 and not verify:
        affected_upload_ids = kwargs.get('affected_upload_ids', [upload_id])
        expected_metadata = kwargs.get('expected_metadata', metadata)

    # Perform edit request
    mer = MetadataEditRequest(
        upload_id=upload_id, query=query, metadata=metadata, entries=entries, verify=verify)
    edit_start = datetime.utcnow().isoformat()[0:22]
    _response, status_code = Upload.edit_metadata(mer, user)
    # Validate result
    assert status_code == expected_status_code, 'Wrong status code returned'
    if status_code == 200 and not verify:
        for upload_id in affected_upload_ids:
            upload = Upload.get(upload_id)
            upload.block_until_complete()
            for entry in upload.calcs:
                assert entry.last_edit_time
                assert entry.last_edit_time.isoformat()[0:22] >= edit_start
                entry_metadata_mongo = entry.mongo_metadata(upload).m_to_dict()
                entry_metadata_es = search(owner=None, query={'calc_id': entry.calc_id}).data[0]
                values_to_check = expected_metadata
                for quantity_name, value_expected in values_to_check.items():
                    # Note, the expected value is provided on the "request format"
                    quantity = _editable_metadata[quantity_name]
                    if quantity_name == 'embargo_length':
                        assert upload.embargo_length == value_expected
                        assert entry_metadata_mongo['embargo_length'] == value_expected
                        assert entry_metadata_es['with_embargo'] == (value_expected > 0)
                    else:
                        value_mongo = entry_metadata_mongo.get(quantity_name)
                        value_es = entry_metadata_es.get(quantity_name)
                        # coauthors and reviewers are not stored in ES. Instead check viewers and writers
                        if quantity_name == 'coauthors':
                            value_es = entry_metadata_es['writers']
                        elif quantity_name == 'reviewers':
                            value_es = entry_metadata_es['viewers']
                        cmp_value_mongo = convert_to_comparable_value(quantity, value_mongo, 'mongo', user)
                        cmp_value_es = convert_to_comparable_value(quantity, value_es, 'es', user)
                        cmp_value_expected = convert_to_comparable_value(quantity, value_expected, 'request', user)
                        # Verify mongo value
                        assert cmp_value_mongo == cmp_value_expected
                        # Verify ES value
                        if quantity_name == 'license':
                            continue  # Not stored indexed by ES
                        elif quantity_name == 'coauthors':
                            # Check that writers == main_author + coauthors
                            assert cmp_value_es == [upload.main_author] + cmp_value_expected
                        elif quantity_name == 'reviewers':
                            # Check that viewers == main_author + coauthors + reviewers
                            assert set(cmp_value_es) == set(
                                [upload.main_author] + (upload.coauthors or []) + cmp_value_expected)
                        else:
                            assert cmp_value_es == cmp_value_expected


def convert_to_comparable_value(quantity, value, from_format, user):
    '''
    Converts `value` from the given source format ('mongo', 'es', 'request')
    to a value that can be compared (user_id for user references, dataset_id
    for datasets, timestamp strings with no more than millisecond precision, etc).
    List quantities are also guaranteed to be converted to lists.
    '''
    if quantity.is_scalar:
        return convert_to_comparable_value_single(quantity, value, from_format, user)
    if type(value) != list:
        value = [value]
    return [convert_to_comparable_value_single(quantity, v, from_format, user) for v in value]


def convert_to_comparable_value_single(quantity, value, format, user):
    if quantity.type in (str, int, float, bool) or isinstance(quantity.type, metainfo.MEnum):
        if value == '' and format == 'request':
            return None
        return value
    elif quantity.type == metainfo.Datetime:
        if not value:
            return None
        return value[0:22]  # Only compare to the millisecond level (mongo's maximal precision)
    elif isinstance(quantity.type, metainfo.Reference):
        # Should be reference
        verify_reference = quantity.type.target_section_def.section_cls
        if verify_reference in [datamodel.User, datamodel.Author]:
            if format == 'mongo':
                return value
            if format == 'es':
                return value['user_id']
            elif format == 'request':
                try:
                    return datamodel.User.get(user_id=value).user_id
                except KeyError:
                    try:
                        return datamodel.User.get(username=value).user_id
                    except KeyError:
                        return datamodel.User.get(email=value).user_id
        elif verify_reference == datamodel.Dataset:
            if format == 'mongo':
                return value
            elif format == 'es':
                return value['dataset_id']
            elif format == 'request':
                try:
                    return datamodel.Dataset.m_def.a_mongo.get(dataset_id=value).dataset_id
                except KeyError:
                    return datamodel.Dataset.m_def.a_mongo.get(
                        user_id=user.user_id, dataset_name=value).dataset_id
    assert False, 'Unhandled type/source'


@pytest.mark.parametrize('kwargs', [
    pytest.param(
        dict(
            metadata=dict(external_db='bad value'),
            expected_status_code=400),
        id='bad-external_db'),
    pytest.param(
        dict(
            metadata=dict(coauthors='silly value'),
            expected_status_code=400),
        id='bad-coauthor-ref'),
    pytest.param(
        dict(
            metadata=dict(reviewers='silly value'),
            expected_status_code=400),
        id='bad-reviewer-ref'),
    pytest.param(
        dict(
            metadata=dict(datasets=['silly value']),
            expected_status_code=400),
        id='bad-dataset-ref'),
    pytest.param(
        dict(
            upload_id='id_published_w',
            metadata=dict(embargo_length=0)),
        id='lift-embargo'),
])
def test_edit_metadata(proc_infra, example_data_writeable, a_dataset, test_users_dict, kwargs):
    kwargs['user'] = test_users_dict[kwargs.get('user', 'test_user')]
    assert_edit_request(**kwargs)


def test_set_and_clear_all(proc_infra, example_data_writeable, a_dataset, test_user):
    # Set all fields a coauthor can set
    assert_edit_request(
        user=test_user,
        metadata=dict(
            upload_name='a humble upload name',
            embargo_length=14,
            coauthors=['lhofstadter'],
            external_id='31415926536',
            comment='a humble comment',
            references=['a reference', 'another reference'],
            external_db='AFLOW',
            reviewers=['lhofstadter'],
            datasets=a_dataset.dataset_id))
    # Clear all fields that can be cleared with a 'set' operation
    # = all of the above, except embargo_length and datasets
    assert_edit_request(
        user=test_user,
        metadata=dict(
            upload_name='',
            coauthors=[],
            external_id='',
            comment='',
            references=[],
            external_db='',
            reviewers=[]))


def test_admin_quantities(proc_infra, example_data_writeable, test_user, other_test_user, admin_user):
    now = datetime.utcnow()
    now_iso = now.isoformat()
    admin_fields = dict(
        upload_create_time=now_iso,
        entry_create_time=now_iso,
        publish_time=now_iso,
        license='a license',
        main_author=other_test_user.user_id)
    assert_edit_request(
        user=admin_user, upload_id='id_published_w', metadata=admin_fields)
    # try to do the same as a non-admin
    for k, v in admin_fields.items():
        assert_edit_request(
            user=test_user, upload_id='id_published_w', metadata={k: v}, expected_status_code=401)
