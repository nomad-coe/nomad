from nomad.config.models.plugins import SchemaPackageEntryPoint
from pydantic import Field


class CountrySchemaPackageEntryPoint(SchemaPackageEntryPoint):
    def load(self):
        from nomad_countries.schema_packages.country import m_package

        return m_package


country = CountrySchemaPackageEntryPoint(
    name='CountryPackage',
    description='Schema package defined using the new plugin mechanism.',
)
