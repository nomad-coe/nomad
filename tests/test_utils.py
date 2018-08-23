import json
import pickle

from nomad.utils import DataObject


class ExampleClass(DataObject):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.msg = 'Hello'


def test_data_objects():
    example = ExampleClass()
    assert example['msg'] == 'Hello'
    assert example.msg == 'Hello'
    json.dumps(example)
    pickled = pickle.dumps(example)
    example = pickle.loads(pickled)
    assert example['msg'] == 'Hello'
    assert example.msg == 'Hello'


def test_data_object_update():
    example = ExampleClass()
    example.update(ExampleClass(more='Hi'))
    assert json.loads(json.dumps(example)).get('more', None) == 'Hi'

    example = ExampleClass()
    example.update({'more': 'Hi'})
    assert json.loads(json.dumps(example)).get('more', None) == 'Hi'

    example.update({'more': None})
    assert example.more == 'Hi'
