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

from fastapi.exceptions import RequestValidationError

from nomad import datamodel, metainfo
from nomad.processing import Upload, MetadataEditRequestHandler
from nomad.processing.data import editable_metadata, mongo_upload_metadata
from nomad.search import search


all_coauthor_metadata = dict(
    # All attributes which a coauthor+ can edit
    upload_name='a humble upload name',
    embargo_length=14,
    coauthors=['lhofstadter'],
    external_id='31415926536',
    comment='a humble comment',
    references=['http://reference1.com', 'http://reference2.com'],
    external_db='AFLOW',
    reviewers=['lhofstadter'],
    datasets=['test_dataset_1'],
)

all_coauthor_upload_metadata = {
    k: v for k, v in all_coauthor_metadata.items() if k in mongo_upload_metadata
}

all_coauthor_entry_metadata = {
    k: v for k, v in all_coauthor_metadata.items() if k not in mongo_upload_metadata
}

all_admin_metadata = dict(
    # Every attribute which only admins can set
    upload_create_time='2021-05-04T00:00:00',
    entry_create_time='2021-05-04T00:00:00',
    publish_time='2021-05-04T00:00:00',
    license='a license',
    main_author='lhofstadter',
)

all_admin_entry_metadata = {
    k: v for k, v in all_admin_metadata.items() if k not in mongo_upload_metadata
}


def assert_edit_request(user, **kwargs):
    # Extract test parameters (lots of defaults)
    upload_id = kwargs.get('upload_id', 'id_unpublished_w')
    query = kwargs.get('query')
    owner = kwargs.get('owner')
    metadata = kwargs.get('metadata')
    entries = kwargs.get('entries')
    entries_key = kwargs.get('entries_key', 'entry_id')
    verify_only = kwargs.get('verify_only', False)
    expected_error_loc = kwargs.get('expected_error_loc')
    affected_upload_ids = kwargs.get('affected_upload_ids', [upload_id])
    expected_metadata = kwargs.get('expected_metadata', metadata)
    # Perform edit request
    edit_request_json = dict(
        query=query,
        owner=owner,
        metadata=metadata,
        entries=entries,
        entries_key=entries_key,
        verify=verify_only,
    )
    edit_start = datetime.utcnow().isoformat()[0:22]
    error_locs = []
    try:
        MetadataEditRequestHandler.edit_metadata(edit_request_json, upload_id, user)
    except RequestValidationError as e:
        error_locs = [error_dict['loc'] for error_dict in e.errors()]
    # Validate result
    if expected_error_loc:
        assert expected_error_loc in error_locs
    if not expected_error_loc and not verify_only:
        assert_metadata_edited(
            user,
            upload_id,
            query,
            metadata,
            entries,
            entries_key,
            verify_only,
            expected_metadata,
            affected_upload_ids,
            edit_start,
        )


def assert_metadata_edited(
    user,
    upload_id,
    query,
    metadata,
    entries,
    entries_key,
    verify_only,
    expected_metadata,
    affected_upload_ids,
    edit_start,
):
    for upload_id in affected_upload_ids:
        upload = Upload.get(upload_id)
        upload.block_until_complete()
        for entry in upload.successful_entries:
            assert entry.last_edit_time
            assert (
                edit_start is None
                or entry.last_edit_time.isoformat()[0:22] >= edit_start
            )
            entry_metadata_mongo = entry.mongo_metadata(upload).m_to_dict()
            entry_metadata_es = search(
                owner=None, query={'entry_id': entry.entry_id}
            ).data[0]
            values_to_check = expected_metadata
            for quantity_name, value_expected in values_to_check.items():
                # Note, the expected value is provided on the "request format"
                quantity = editable_metadata[quantity_name]
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
                    cmp_value_mongo = convert_to_comparable_value(
                        quantity, value_mongo, 'mongo', user
                    )
                    cmp_value_es = convert_to_comparable_value(
                        quantity, value_es, 'es', user
                    )
                    cmp_value_expected = convert_to_comparable_value(
                        quantity, value_expected, 'request', user
                    )
                    # Verify mongo value
                    assert (
                        cmp_value_mongo == cmp_value_expected
                    ), f'Wrong mongo value for {quantity_name}'
                    # Verify ES value
                    if quantity_name == 'license':
                        continue  # Not stored indexed by ES
                    elif quantity_name == 'coauthors':
                        # Check that writers == main_author + coauthors
                        assert (
                            cmp_value_es == [upload.main_author] + cmp_value_expected
                        ), f'Wrong es value for {quantity_name}'
                    elif quantity_name == 'reviewers':
                        # Check that viewers == main_author + coauthors + reviewers
                        assert set(cmp_value_es) == set(
                            [upload.main_author]
                            + (upload.coauthors or [])
                            + cmp_value_expected
                        ), f'Wrong es value for {quantity_name}'
                    else:
                        assert (
                            cmp_value_es == cmp_value_expected
                        ), f'Wrong es value for {quantity_name}'


def convert_to_comparable_value(quantity, value, from_format, user):
    """
    Converts `value` from the given source format ('mongo', 'es', 'request')
    to a value that can be compared (user_id for user references, dataset_id
    for datasets, timestamp strings with no more than millisecond precision, etc).
    List quantities are also guaranteed to be converted to lists.
    """
    if quantity.is_scalar:
        return convert_to_comparable_value_single(quantity, value, from_format, user)
    if value is None and from_format == 'es':
        return []
    if not isinstance(value, list):
        value = [value]
    return [
        convert_to_comparable_value_single(quantity, v, from_format, user)
        for v in value
    ]


def convert_to_comparable_value_single(quantity, value, format, user):
    if quantity.type in (str, int, float, bool) or isinstance(
        quantity.type, metainfo.MEnum
    ):
        if value == '' and format == 'request':
            return None
        return value
    elif quantity.type == metainfo.Datetime:
        if not value:
            return None
        # datetime that is returned from Datetime class now has the timezone information.
        # Hence, need to drop the last 3 elements in the datetime string [:19] instead of [:22]
        return value[
            0:19
        ]  # Only compare to the millisecond level (mongo's maximal precision).
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
                    return datamodel.Dataset.m_def.a_mongo.get(
                        dataset_id=value
                    ).dataset_id
                except KeyError:
                    return datamodel.Dataset.m_def.a_mongo.get(
                        user_id=user.user_id, dataset_name=value
                    ).dataset_id
    assert False, 'Unhandled type/source'


@pytest.mark.parametrize(
    'kwargs',
    [
        pytest.param(
            dict(
                metadata=dict(external_db='bad value'),
                expected_error_loc=('metadata', 'external_db'),
            ),
            id='bad-external_db',
        ),
        pytest.param(
            dict(
                metadata=dict(coauthors='silly value'),
                expected_error_loc=('metadata', 'coauthors'),
            ),
            id='bad-coauthor-ref',
        ),
        pytest.param(
            dict(
                metadata=dict(reviewers='silly value'),
                expected_error_loc=('metadata', 'reviewers'),
            ),
            id='bad-reviewer-ref',
        ),
        pytest.param(
            dict(
                metadata=dict(datasets=['silly value']),
                expected_error_loc=('metadata', 'datasets'),
            ),
            id='bad-dataset-ref',
        ),
        pytest.param(
            dict(upload_id='id_published_w', metadata=dict(embargo_length=0)),
            id='lift-embargo',
        ),
        pytest.param(
            dict(
                query={
                    'and': [
                        {'upload_create_time:gt': '2021-01-01'},
                        {'published': False},
                    ]
                },
                owner='user',
                upload_id=None,
                metadata=dict(comment='new comment'),
                affected_upload_ids=['id_unpublished_w'],
            ),
            id='query-ok',
        ),
        pytest.param(
            dict(
                query={'upload_create_time:lt': '2021-01-01'},
                owner='user',
                upload_id=None,
                metadata=dict(comment='new comment'),
                expected_error_loc=('query',),
            ),
            id='query-no-results',
        ),
        pytest.param(
            dict(
                query={'upload_create_time:gt': '2021-01-01'},
                owner='user',
                upload_id=None,
                metadata=dict(comment='new comment'),
                affected_upload_ids=['id_unpublished_w'],
            ),
            id='query-contains-published',
        ),
    ],
)
def test_edit_metadata(
    proc_infra,
    purged_app,
    example_data_writeable,
    example_datasets,
    test_users_dict,
    kwargs,
):
    kwargs['user'] = test_users_dict[kwargs.get('user', 'test_user')]
    assert_edit_request(**kwargs)


def test_set_and_clear_all(
    proc_infra, example_data_writeable, example_datasets, test_user
):
    # Set all fields a coauthor can set
    assert_edit_request(user=test_user, metadata=all_coauthor_metadata)
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
            reviewers=[],
        ),
    )


@pytest.mark.parametrize(
    'kwargs',
    [
        pytest.param(
            dict(
                metadata_1=dict(datasets={'add': 'test_dataset_1'}),
                expected_metadata_1=dict(datasets=['test_dataset_1']),
                metadata_2=dict(datasets={'add': 'test_dataset_2'}),
                expected_metadata_2=dict(datasets=['test_dataset_1', 'test_dataset_2']),
            ),
            id='datasets-add',
        ),
        pytest.param(
            dict(
                metadata_1=dict(datasets=['test_dataset_1']),
                metadata_2=dict(datasets=['ref:test_dataset_2']),
                expected_metadata_2=dict(datasets=['test_dataset_1', 'test_dataset_2']),
            ),
            id='datasets-implicit-add',
        ),
        pytest.param(
            dict(
                metadata_1=dict(datasets=['test_dataset_1']),
                metadata_2=dict(datasets=['test_dataset_3']),
                expected_error_loc_2=('metadata', 'datasets'),
            ),
            id='datasets-implicit-add-not-owner',
        ),
        pytest.param(
            dict(
                metadata_1=dict(datasets=['test_dataset_1']),
                metadata_2=dict(datasets=['ref:test_dataset_3']),
                expected_error_loc_2=('metadata', 'datasets'),
            ),
            id='datasets-implicit-add-not-owner-ref',
        ),
        pytest.param(
            dict(
                metadata_1=dict(datasets=['test_dataset_1', 'test_dataset_2']),
                metadata_2=dict(datasets={'remove': 'test_dataset_1'}),
                expected_metadata_2=dict(datasets=['test_dataset_2']),
            ),
            id='datasets-remove',
        ),
        pytest.param(
            dict(
                metadata_1=dict(datasets=['test_dataset_1', 'ref:test_dataset_2']),
                metadata_2=dict(datasets={'remove': 'test_dataset_2'}),
                expected_error_loc_2=('metadata', 'datasets'),
            ),
            id='datasets-remove-with-doi',
        ),
        pytest.param(
            dict(
                metadata_1=dict(datasets='test_dataset_1'),
                metadata_2=dict(
                    datasets={'add': ['test_dataset_2'], 'remove': 'ref:test_dataset_1'}
                ),
                expected_metadata_2=dict(datasets=['test_dataset_2']),
            ),
            id='datasets-add+remove',
        ),
        pytest.param(
            dict(
                metadata_1=dict(datasets={'set': ['test_dataset_2']}),
                expected_error_loc_1=('metadata', 'datasets'),
                metadata_2=dict(
                    datasets={'add': 'test_dataset_1', 'set': ['test_dataset_2']}
                ),
                expected_error_loc_2=('metadata', 'datasets'),
            ),
            id='datasets-set',
        ),
        pytest.param(
            dict(
                metadata_1=dict(coauthors='scooper'),
                expected_metadata_1=dict(coauthors=[]),
                metadata_2=dict(coauthors={'add': ['lhofstadter']}),
                expected_metadata_2=dict(coauthors=['lhofstadter']),
            ),
            id='coauthors-add',
        ),
        pytest.param(
            dict(
                metadata_1=dict(coauthors='lhofstadter'),
                expected_metadata_1=dict(coauthors=['lhofstadter']),
                metadata_2=dict(coauthors={'add': ['admin'], 'remove': 'lhofstadter'}),
                expected_metadata_2=dict(coauthors=['admin']),
            ),
            id='coauthors-add+remove',
        ),
        pytest.param(
            dict(
                metadata_1=dict(coauthors='lhofstadter'),
                expected_metadata_1=dict(coauthors=['lhofstadter']),
                metadata_2=dict(coauthors={'add': ['admin'], 'set': 'lhofstadter'}),
                expected_error_loc_2=('metadata', 'coauthors'),
            ),
            id='coauthors-add+set',
        ),
        pytest.param(
            dict(
                metadata_1=dict(reviewers='scooper'),
                expected_metadata_1=dict(reviewers=[]),
                metadata_2=dict(reviewers={'add': ['lhofstadter']}),
                expected_metadata_2=dict(reviewers=['lhofstadter']),
            ),
            id='reviewers-add',
        ),
        pytest.param(
            dict(
                metadata_1=dict(reviewers='lhofstadter'),
                expected_metadata_1=dict(reviewers=['lhofstadter']),
                metadata_2=dict(reviewers={'add': ['admin'], 'remove': 'lhofstadter'}),
                expected_metadata_2=dict(reviewers=['admin']),
            ),
            id='reviewers-add+remove',
        ),
        pytest.param(
            dict(
                metadata_1=dict(reviewers='lhofstadter'),
                expected_metadata_1=dict(reviewers=['lhofstadter']),
                metadata_2=dict(reviewers={'remove': ['admin'], 'set': 'lhofstadter'}),
                expected_error_loc_2=('metadata', 'reviewers'),
            ),
            id='reviewers-remove+set',
        ),
        pytest.param(
            dict(
                metadata_1=dict(references='http://ref1.com'),
                expected_metadata_1=dict(references=['http://ref1.com']),
                metadata_2=dict(references={'add': ['http://ref2.com']}),
                expected_metadata_2=dict(
                    references=['http://ref1.com', 'http://ref2.com']
                ),
            ),
            id='references-add',
        ),
        pytest.param(
            dict(
                metadata_1=dict(references=['http://ref1.com', 'http://ref2.com']),
                expected_metadata_1=dict(
                    references=['http://ref1.com', 'http://ref2.com']
                ),
                metadata_2=dict(
                    references={'add': 'http://ref3.com', 'remove': ['http://ref1.com']}
                ),
                expected_metadata_2=dict(
                    references=['http://ref2.com', 'http://ref3.com']
                ),
            ),
            id='references-add+remove',
        ),
        pytest.param(
            dict(
                metadata_1=dict(references='http://ref1.com'),
                expected_metadata_1=dict(references=['http://ref1.com']),
                metadata_2=dict(
                    references={
                        'remove': 'http://ref4.com',
                        'add': [
                            'http://ref2.com',
                            'http://ref3.com',
                            'http://ref4.com',
                        ],
                    }
                ),
                expected_error_loc_2=('metadata', 'references'),
            ),
            id='references-add+remove-incoherent',
        ),
        pytest.param(
            dict(
                metadata_1=dict(references='http://ref1.com'),
                expected_metadata_1=dict(references=['http://ref1.com']),
                metadata_2=dict(references={'add': ['http://ref2', 'http://ref3.com']}),
                expected_error_loc_2=('metadata', 'references'),
            ),
            id='references-not-valid-URL',
        ),
    ],
)
def test_list_quantities(
    proc_infra,
    purged_app,
    example_data_writeable,
    example_datasets,
    test_users_dict,
    kwargs,
):
    def replace_dataset_ref(dataset_ref):
        if dataset_ref.startswith('ref:'):
            for dataset in example_datasets:
                if dataset.dataset_name == dataset_ref[4:]:
                    return dataset.dataset_id
        return dataset_ref

    def replace_dataset_ref_or_reflist(ref_or_reflist):
        if isinstance(ref_or_reflist, list):
            return [replace_dataset_ref(ref) for ref in ref_or_reflist]
        return replace_dataset_ref(ref_or_reflist)

    kwargs['user'] = test_users_dict[kwargs.get('user', 'test_user')]
    for suffix in ('_1', '_2'):
        for arg in ('metadata', 'expected_metadata', 'expected_error_loc'):
            if arg in kwargs:
                kwargs.pop(arg)
            if arg + suffix in kwargs:
                kwargs[arg] = kwargs.pop(arg + suffix)
        datasets = kwargs['metadata'].get('datasets')
        if datasets is not None:
            if isinstance(datasets, dict):
                datasets = {
                    op: replace_dataset_ref_or_reflist(v) for op, v in datasets.items()
                }
            else:
                datasets = replace_dataset_ref_or_reflist(datasets)
            kwargs['metadata']['datasets'] = datasets
        assert_edit_request(**kwargs)


def test_admin_quantities(
    proc_infra, example_data_writeable, test_user, other_test_user, admin_user
):
    assert_edit_request(
        user=admin_user, upload_id='id_published_w', metadata=all_admin_metadata
    )
    # try to do the same as a non-admin
    for k, v in all_admin_metadata.items():
        assert_edit_request(
            user=test_user,
            upload_id='id_unpublished_w',
            metadata={k: v},
            expected_error_loc=('metadata', k),
        )


def test_query_cannot_set_upload_attributes(
    proc_infra, example_data_writeable, example_datasets, test_user
):
    query = {'and': [{'upload_create_time:gt': '2021-01-01'}, {'published': False}]}
    for k, v in all_coauthor_upload_metadata.items():
        # Attempting to edit an upload level attribute with query should always fail,
        # regardless of if upload_id is specified
        for upload_id in (None, 'id_unpublished_w'):
            assert_edit_request(
                user=test_user,
                query=query,
                owner='user',
                upload_id=upload_id,
                metadata={k: v},
                expected_error_loc=('metadata', k),
            )
    # Attempting to edit an entry level attribute with query should always succeed
    assert_edit_request(
        user=test_user,
        query=query,
        owner='user',
        upload_id=None,
        metadata=all_coauthor_entry_metadata,
        affected_upload_ids=['id_unpublished_w'],
    )
