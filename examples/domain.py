from nomad import datamodel

print(datamodel.EntryMetadata(domain='DFT', calc_id='test').__class__.__name__)
print(datamodel.EntryMetadata(calc_id='test').__class__.__name__)
print(datamodel.EntryMetadata(domain='EMS', calc_id='test').__class__.__name__)
