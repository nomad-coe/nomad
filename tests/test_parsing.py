from nomad.parsing import JSONStreamWriter
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
    writer = JSONStreamWriter(out, pretty=pretty)
    writer.open_array()
    writer.open_object()
    writer.key('key1')
    writer.value('value')
    writer.key('key2')
    writer.value(1)
    writer.close_object()
    writer.open_object()
    writer.key('key')
    writer.open_object()
    writer.key('key')
    writer.value('value')
    writer.close_object()
    writer.close_object()
    writer.close_array()
    writer.close()

    assert create_reference(example_data, pretty) == out.getvalue()

