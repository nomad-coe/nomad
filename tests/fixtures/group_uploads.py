"""
Fields:
- C: coauthors (users)
- R: reviewers (users)
- CG: coauthor_groups (agents)
- RG: reviewer_groups (agents)

Values:
- XMN: field X with userM and userN
- XgKLgMN: field X with groupKL and groupMN
"""

from typing import Sequence
import pytest

from nomad.utils.exampledata import ExampleData


@pytest.fixture('session')
def group_upload_molds(
    convert_user_labels_to_ids, convert_group_labels_to_ids, user1, user2
):
    default = {'main_author': user1}
    molds = {
        'no_group': {},
        'user2': {'main_author': user2},
        'embargo0': {'published': True, 'embargo_length': 0},
        'embargo3': {'published': True, 'embargo_length': 3},
        'C2': {'coauthors': ['user2']},
        'R2': {'reviewers': ['user2']},
        'CGg1': {'coauthor_groups': ['group1']},
        'CGg2': {'coauthor_groups': ['group2']},
        'CGg123': {'coauthor_groups': ['group123']},
        'RGg1': {'reviewer_groups': ['group1']},
        'RGg2': {'reviewer_groups': ['group2']},
        'RGg123': {'reviewer_groups': ['group123']},
        'RGall': {'reviewer_groups': ['all']},
        'full_agents': {
            'coauthors': ['user2', 'user4'],
            'reviewers': ['user3', 'user5', 'user7'],
            'coauthor_groups': ['group2', 'group8', 'group14'],
            'reviewer_groups': ['group3', 'group8', 'group15'],
        },
    }
    molds = {k: {**default, **v} for k, v in molds.items()}

    for label, mold in molds.items():
        mold['upload_id'] = f'id_{label}'

        convert = convert_user_labels_to_ids
        mold['coauthors'] = convert(mold.get('coauthors', []))
        mold['reviewers'] = convert(mold.get('reviewers', []))

        convert = convert_group_labels_to_ids
        mold['coauthor_groups'] = convert(mold.get('coauthor_groups', []))
        mold['reviewer_groups'] = convert(mold.get('reviewer_groups', []))

    return molds


@pytest.fixture('session')
def create_group_uploads_from_molds(group_upload_molds):
    """Returned function creates and returns uploads with given labels.

    Use `yield from create_group_uploads_from_molds(...)` for the actual upload fixtures.
    """

    def create(labels: Sequence[str], save_kwargs=None):
        data = ExampleData()
        for label in labels:
            mold = group_upload_molds[label]
            upload_id = mold['upload_id']
            entry_id = f'{upload_id}_1'

            data.create_upload(**mold)
            data.create_entry(upload_id=upload_id, entry_id=entry_id)

        if save_kwargs is None:
            save_kwargs = {}

        data.save(**save_kwargs)
        yield data
        data.delete()

    return create


@pytest.fixture(scope='module')
def uploads_get_groups(create_group_uploads_from_molds, elastic_module, groups_module):
    labels = ('no_group', 'CGg2', 'RGg2', 'CGg123', 'RGg123', 'RGall', 'CG2', 'RG2')
    yield from create_group_uploads_from_molds(labels)


@pytest.fixture(scope='function')
def upload_no_group(create_group_uploads_from_molds, groups_function):
    labels = ('no_group',)
    yield from create_group_uploads_from_molds(labels)


@pytest.fixture(scope='function')
def upload_CGg2g123(
    create_group_uploads_from_molds,
    groups_function,
):
    labels = ('CGg2g123',)
    yield from create_group_uploads_from_molds(labels)


@pytest.fixture(scope='function')
def upload_RGall(create_group_uploads_from_molds, mongo_function):
    labels = ('RGall',)
    yield from create_group_uploads_from_molds(labels)


@pytest.fixture(scope='class')
def uploads_search_query_groups(
    create_group_uploads_from_molds,
    elastic_module,
    groups_module,
):
    labels = (
        'embargo0',
        'embargo3',
        'no_group',
        'CGg1',
        'RGg1',
        'CGg2',
        'RGg2',
        'CGg123',
        'RGg123',
        'RGall',
        'user2',
    )
    yield from create_group_uploads_from_molds(
        labels, save_kwargs={'with_files': False, 'with_mongo': False}
    )
