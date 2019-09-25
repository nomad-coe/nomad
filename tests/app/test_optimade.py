import json

from nomad.processing import Upload
from nomad.search import SearchRequest


def test_get_entry(published: Upload):
    calc_id = list(published.calcs)[0].calc_id

    with published.upload_files.archive_file(calc_id) as f:
        data = json.load(f)

    assert 'OptimadeStructureEntry' in data
    search_result = SearchRequest().search_parameter('calc_id', calc_id).execute_paginated()['results'][0]
    assert 'optimade' in search_result
