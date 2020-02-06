from nomad import datamodel

print(datamodel.CalcWithMetadata(domain='DFT', calc_id='test').__class__.__name__)
print(datamodel.CalcWithMetadata(calc_id='test').__class__.__name__)
print(datamodel.CalcWithMetadata(domain='EMS', calc_id='test').__class__.__name__)
