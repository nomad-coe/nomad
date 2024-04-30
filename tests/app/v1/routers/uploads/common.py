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


def assert_entry(entry, has_metadata=True, **kwargs):
    """Checks the content of a returned entry dictionary."""
    assert 'upload_id' in entry
    assert 'entry_id' in entry
    assert 'entry_create_time' in entry
    assert not entry['process_running']
    for key, value in kwargs.items():
        assert entry.get(key, None) == value
    if has_metadata:
        assert 'entry_metadata' in entry
