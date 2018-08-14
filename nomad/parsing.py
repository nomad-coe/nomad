import json


class JSONStreamWriter():
    START = 0
    OBJECT = 1
    ARRAY = 2
    KEY_VALUE = 3

    """
    A generator that allows to output JSON based on calling 'event' functions.
    Its pure python and could be replaced by some faster implementation, e.g. yajl-py.
    It uses standard json decode to write values. This allows to mix streaming with
    normal encoding.

    Arguments:
        file: A file like to write to.
        pretty: True to indent and use separators.

    Raises:
        AssertionError: If methods were called in a non JSON fashion. Call :func:`close`
        to make sure everything was closed properly.
    """
    def __init__(self, file, pretty=False):
        self._fp = file
        self._pretty = pretty

        self._indent = ''  # the current indent
        self._separators = ['']  # a stack of the next necessary separator
        self._states = [JSONStreamWriter.START]  # a stack of what is currenty open

    def _write(self, str):
        self._fp.write(str)

    def _write_seperator(self):
        self._write(self._separators.pop())

    def _seperator_with_newline(self, base=None):
        pretty_ext = ('\n%s' % self._indent) if self._pretty else ''
        if base is None:
            return pretty_ext
        else:
            return '%s%s' % (base, pretty_ext)

    def _open(self, open_char):
        self._write_seperator()
        self._write(open_char)
        self._indent = '%s  ' % self._indent
        self._separators.append(self._seperator_with_newline())

    def _close(self, close_char):
        self._separators.pop()
        self._indent = self._indent[:-2]
        self._write(self._seperator_with_newline())
        self._write(close_char)
        self._separators.append(self._seperator_with_newline(','))

    def open_object(self):
        assert self._states[-1] != JSONStreamWriter.OBJECT, "Cannot open object in object."
        if self._states[-1] == JSONStreamWriter.KEY_VALUE:
            self._states.pop()
        self._open('{')
        self._states.append(JSONStreamWriter.OBJECT)

    def close_object(self):
        assert self._states.pop() == JSONStreamWriter.OBJECT, "Can only close object in object."
        self._close('}')

    def open_array(self):
        assert self._states[-1] != JSONStreamWriter.OBJECT, "Cannot open array in object."
        if self._states[-1] == JSONStreamWriter.KEY_VALUE:
            self._states.pop()
        self._open('[')
        self._states.append(JSONStreamWriter.ARRAY)

    def close_array(self):
        assert self._states.pop() == JSONStreamWriter.ARRAY, "Can only close array in array."
        self._close(']')

    def key_value(self, key, value):
        self.key(key)
        self.value(value)

    def key(self, key):
        assert self._states[-1] == JSONStreamWriter.OBJECT, "Key can only be in objects."
        self._write_seperator()
        json.dump(key, self._fp)
        self._separators.append(': ' if self._pretty else ':')
        self._states.append(JSONStreamWriter.KEY_VALUE)

    def value(self, value):
        assert self._states[-1] != JSONStreamWriter.OBJECT, "Values can not be in objects."
        if self._states[-1] == JSONStreamWriter.KEY_VALUE:
            self._states.pop()

        self._write_seperator()
        json.dump(value, self._fp)
        self._separators.append(self._seperator_with_newline(','))

    def close(self):
        assert self._states[-1] == JSONStreamWriter.START, "Something was not closed."
