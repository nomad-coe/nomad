import json


class JSONStreamGenerator():
    """
    A generator that allows to output JSON based on calling 'event' functions.
    Its pure python and could be replaced by some faster implementation, e.g. yajl-py.
    It does not do anychecks. Expect random exceptions when events are out of order or
    imcomplete.
    """
    def __init__(self, fp, pretty=False):
        self._fp = fp
        self._pretty = pretty

        self._indent = ''
        self._separators = ['']

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
        self._open('{')

    def close_object(self):
        self._close('}')

    def open_array(self):
        self._open('[')

    def close_array(self):
        self._close(']')

    def key_value(self, key, value):
        self.key(key)
        self.value(value)

    def key(self, key):
        self._write_seperator()
        json.dump(key, self._fp)
        self._separators.append(': ' if self._pretty else ':')

    def value(self, value):
        self._write_seperator()
        json.dump(value, self._fp)
        self._separators.append(self._seperator_with_newline(','))
