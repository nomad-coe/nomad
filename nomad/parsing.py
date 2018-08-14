from abc import ABCMeta, abstractmethod
import json

from nomadcore.local_backend import LocalBackend as LegacyLocalBackend
from nomadcore.local_backend import Section


class DelegatingMeta(ABCMeta):
    def __new__(meta, name, bases, dct):
        abstract_method_names = frozenset.union(*(base.__abstractmethods__ for base in bases))
        for name in abstract_method_names:
            if name not in dct:
                dct[name] = DelegatingMeta._make_delegator_method(name)

        return super(DelegatingMeta, meta).__new__(meta, name, bases, dct)

    @staticmethod
    def _make_delegator_method(name):
        def delegator(self, *args, **kwargs):
            return getattr(self._delegate, name)(*args, **kwargs)
        return delegator


class AbstractParserBackend(metaclass=ABCMeta):

    @abstractmethod
    def metaInfoEnv(self):
        """ Returns the meta info used by this backend. """
        pass

    @abstractmethod
    def startedParsingSession(
            self, mainFileUri, parserInfo, parserStatus=None, parserErrors=None):
        """
        Should be called when the parsing starts.
        ParserInfo should be a valid json dictionary.
        """
        pass

    @abstractmethod
    def finishedParsingSession(
            self, parserStatus, parserErrors, mainFileUri=None, parserInfo=None,
            parsingStats=None):
        """ Called when the parsing finishes. """
        pass

    @abstractmethod
    def openSection(self, metaName):
        """ Opens a new section and returns its new unique gIndex. """
        pass

    @abstractmethod
    def closeSection(self, metaName, gIndex):
        """
        Closes the section with the given meta name and index. After this, no more
        value can be added to this section.
        """
        pass

    @abstractmethod
    def openNonOverlappingSection(self, metaName):
        """ Opens a new non overlapping section. """
        pass

    @abstractmethod
    def closeNonOverlappingSection(self, metaName):
        """
        Closes the current non overlapping section for the given meta name. After
        this, no more value can be added to this section.
        """
        pass

    @abstractmethod
    def openSections(self):
        """ Returns the sections that are still open as metaName, gIndex tuples. """
        pass

    @abstractmethod
    def addValue(self, metaName, value, gIndex=-1):
        """
        Adds a json value for the given metaName. The gIndex is used to identify
        the right parent section.
        """
        pass

    @abstractmethod
    def addRealValue(self, metaName, value, gIndex=-1):
        """
        Adds a float value for the given metaName. The gIndex is used to identify
        the right parent section.
        """
        pass

    @abstractmethod
    def addArray(self, metaName, shape, gIndex=-1):
        """
        Adds an unannitialized array of the given shape for the given metaName.
        The gIndex is used to identify the right parent section.
        This is neccessary before array values can be set with :func:`setArrayValues`.
        """

    @abstractmethod
    def setArrayValues(self, metaName, values, offset=None, gIndex=-1):
        """
        Adds values of the given numpy array to the last array added for the given
        metaName and parent gIndex.
        """
        pass

    @abstractmethod
    def addArrayValues(self, metaName, values, gIndex=-1):
        """
        Adds an array with the given numpy array values for the given metaName and
        parent section gIndex.
        """
        pass


class LegacyParserBackend(AbstractParserBackend, metaclass=DelegatingMeta):
    """
    Simple implementation of :class:`AbstractParserBackend` that delegates all calls to
    another parser object that not necessarely need to decend from the abstract base class.
    """
    def __init__(self, legacy_backend):
        self._delegate = legacy_backend


class LocalBackend(LegacyParserBackend):
    def __init__(self, *args, **kwargs):
        delegate = LegacyLocalBackend(*args, **kwargs)
        super().__init__(delegate)

    @staticmethod
    def _write(json_writer, value):
        if isinstance(value, list):
            json_writer.open_array()
            for item in value:
                LocalBackend._write(json_writer, item)
            json_writer.close_array()
        elif isinstance(value, Section):
            section = value
            json_writer.open_object()
            json_writer.key_value('_name', section.name)
            json_writer.key_value('_gIndex', section.gIndex)
            for name, value in section.items():
                json_writer.key(name)
                LocalBackend._write(json_writer, value)
            json_writer.close_object()
        else:
            json_writer.value(value)

    def write_json(self, json_writer):
        json_writer.open_object()
        json_writer.key('section_run')
        json_writer.open_array()
        for run in self._delegate.results['section_run']:
            LocalBackend._write(json_writer, run)
        json_writer.close_array()
        json_writer.close_object()
        json_writer.close()

        self._delegate.results.print_summary()


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
