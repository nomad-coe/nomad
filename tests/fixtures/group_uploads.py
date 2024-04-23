import pytest

from nomad.utils.exampledata import ExampleData


@pytest.fixture('module')
def fill_group_data(convert_group_labels_to_ids):
    def fill(data: ExampleData, id_label, c_groups, r_groups, **kwargs):
        upload_id = f'id_{id_label}'
        entry_id = f'{upload_id}_1'

        d = kwargs.copy()
        d['upload_id'] = upload_id
        d['coauthor_groups'] = convert_group_labels_to_ids(c_groups)
        d['reviewer_groups'] = convert_group_labels_to_ids(r_groups)
        d.setdefault('published', d.get('embargo_length') is not None)
        if d.get('embargo_length') is None:
            d['embargo_length'] = 0

        data.create_upload(**d)
        data.create_entry(upload_id=upload_id, entry_id=entry_id)

    return fill


@pytest.fixture(scope='module')
def uploads_get_groups(
    elastic_module, mongo_module, user_groups_module, user1, fill_group_data
):
    data = ExampleData(main_author=user1)
    fill_group_data(data, 'no_group', [], [])
    fill_group_data(data, 'coauthor_group2', ['group2'], [])
    fill_group_data(data, 'reviewer_group2', [], ['group2'])
    fill_group_data(data, 'coauthor_group012', ['group012'], [])
    fill_group_data(data, 'reviewer_group012', [], ['group012'])
    fill_group_data(data, 'reviewer_all', [], ['all'])

    data.save(with_files=False)

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_no_group(mongo_function, user1):
    data = ExampleData(main_author=user1)
    data.create_upload(upload_id='id_no_group')
    data.save()

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_coauthor_group2_and_group012(
    mongo_function,
    user1,
    group2,
    group012,
):
    data = ExampleData(main_author=user1)
    data.create_upload(
        upload_id='id_coauthor_ogroup_mgroup',
        coauthor_groups=[group2.group_id, group012.group_id],
    )
    data.save()

    yield data

    data.delete()


@pytest.fixture(scope='function')
def upload_reviewer_all_group(mongo_function, user1):
    data = ExampleData(main_author=user1)
    data.create_upload(upload_id='id_reviewer_all', reviewer_groups=['all'])
    data.save()

    yield data

    data.delete()


@pytest.fixture(scope='class')
def uploads_search_query_groups(
    elastic_module, user_groups_module, user1, user2, fill_group_data
):
    data = ExampleData(main_author=user1)
    fill_group_data(data, 'no_embargo', [], [], embargo_length=0)
    fill_group_data(data, 'with_embargo', [], [], embargo_length=3)
    fill_group_data(data, 'no_group', [], [])
    fill_group_data(data, 'coauthor_group1', ['group1'], [])
    fill_group_data(data, 'reviewer_group1', [], ['group1'])
    fill_group_data(data, 'coauthor_group2', ['group2'], [])
    fill_group_data(data, 'reviewer_group2', [], ['group2'])
    fill_group_data(data, 'coauthor_group012', ['group012'], [])
    fill_group_data(data, 'reviewer_group012', [], ['group012'])
    fill_group_data(data, 'reviewer_all', [], ['all'])
    fill_group_data(data, 'user2', [], [], main_author=user2)
    data.save(with_files=False, with_mongo=False)

    yield data

    data.delete()
