def assert_upload(response_json, **kwargs):
    data = response_json['data']
    assert 'upload_id' in response_json
    assert 'upload_id' in data
    assert 'upload_create_time' in data
    assert 'main_author' in data
    assert 'coauthors' in data
    assert 'reviewers' in data
    assert 'viewers' in data
    assert 'writers' in data
    assert 'published' in data
    assert 'with_embargo' in data
    assert 'embargo_length' in data
    assert 'license' in data
    assert (data['embargo_length'] > 0) == data['with_embargo']
    if data['published']:
        assert 'publish_time' in data

    for key, value in kwargs.items():
        assert data.get(key, None) == value
    return data