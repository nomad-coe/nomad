from nomad.parsing import JSONStreamGenerator
from io import StringIO
import json
import pytest


def create_reference(data, pretty):
    if (pretty):
        return json.dumps(data, indent=2)
    else:
        return json.dumps(data, separators=(',', ':'))


@pytest.mark.parametrize("pretty", [False, True])
def test_stream_generator(pretty):
    example_data = [
        {
            'key1': 'value',
            'key2': 1
        },
        {
            'key': {
                'key': 'value'
            }
        }
    ]

    out = StringIO()
    generator = JSONStreamGenerator(out, pretty=pretty)
    generator.open_array()
    generator.open_object()
    generator.key('key1')
    generator.value('value')
    generator.key('key2')
    generator.value(1)
    generator.close_object()
    generator.open_object()
    generator.key('key')
    generator.open_object()
    generator.key('key')
    generator.value('value')
    generator.close_object()
    generator.close_object()
    generator.close_array()

    assert create_reference(example_data, pretty) == out.getvalue()

