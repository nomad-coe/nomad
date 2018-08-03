import abc
import re

from nomadcore.parser_backend import JsonParseEventsWriterBackend
from vaspparser import VASPParser

_meta_info_path = './submodules/nomad-meta-info/meta_info/nomad_meta_info/'


class Parser(abc.ABC):
    """
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.
    """
    def __init__(self, name, main_file_re, main_contents_re):
        super().__init__()
        self.name = name
        self._main_file_re = re.compile(main_file_re)
        self._main_contents_re = re.compile(main_contents_re)

    def is_mainfile(self, upload, filename):
        if self._main_file_re.match(filename):
            file = None
            try:
                file = upload.open_file(filename)
                return self._main_contents_re.match(file.read(500))
            finally:
                if file:
                    file.close()

    @abc.abstractmethod
    def run(self, mainfile):
        pass


class VASPRunParser(Parser):
    def __init__(self):
        super().__init__(
            name='VASPRunParser',
            main_file_re=r'^.*\.xml$',
            main_contents_re=(
                r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
                r'?\s*<modeling>'
                r'?\s*<generator>'
                r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
                r'?')
        )

    def run(self, mainfile):
        parser = VASPParser(backend=JsonParseEventsWriterBackend)
        parser.parse(mainfile)

parsers = [
    VASPRunParser()
]
parser_dict = {parser.name: parser for parser in parsers}
