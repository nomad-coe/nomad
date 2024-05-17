from nomad.config.models.plugins import ParserEntryPoint
from pydantic import Field


class CountryParserEntryPoint(ParserEntryPoint):
    def load(self):
        from nomad_countries.parsers.country import CountryParser

        return CountryParser(**self.dict())


country = CountryParserEntryPoint(
    name='CountryParser',
    description='Parser defined using the new plugin mechanism.',
    mainfile_name_re='^.*\.data\.txt$',
    mainfile_contents_re='#Country=\w+',
)
