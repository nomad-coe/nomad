import json

from nomad import config


def test_mapping_compatibility(elastic_infra):
    from nomad.infrastructure import elastic_client

    v0 = elastic_client.indices.get(config.elastic.index_name)
    v1 = elastic_client.indices.get(config.elastic.entries_index)

    def get_mapping(index):
        assert len(index) == 1
        index = index[next(iter(index))]
        assert len(index['mappings']) == 1
        return index['mappings'][next(iter(index['mappings']))]

    v0, v1 = get_mapping(v0), get_mapping(v1)

    def compare(a, b, path='', results=None):
        if results is None:
            results = []
        if path != '':
            path += '.'
        for key in set(list(a.keys()) + list(b.keys())):
            if key in a and key in b:
                next_a, next_b = a[key], b[key]
                if isinstance(next_a, dict) and isinstance(next_b, dict):
                    compare(next_a, next_b, f'{path}{key}', results=results)
                    continue

                if next_a == next_b:
                    continue

            results.append(f"{'v0' if key in a else 'v1'}:{path}{key}")

        return results

    for diff in compare(v0, v1):
        # assert that there are only top-level differences and mapping types and fields are
        # the same
        assert len([c for c in diff if c == '.']) == 1, diff
