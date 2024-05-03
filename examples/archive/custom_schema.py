from nomad.datamodel import Schema, ArchiveSection
from nomad.metainfo.metainfo import Quantity, Datetime, SubSection


class Sample(ArchiveSection):
    added_date = Quantity(type=Datetime)
    formula = Quantity(type=str)

    sample_id = Quantity(type=str)

    def normalize(self, archive, logger):
        super(Sample, self).normalize(archive, logger)

        if self.sample_id is None:
            self.sample_id = f'{self.added_date}--{self.formula}'


class SampleDatabase(Schema):
    samples = SubSection(section=Sample, repeats=True)
