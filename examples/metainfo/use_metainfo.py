import json

from nomad.datamodel.metainfo import public

# A simple example that demonstrates how to set references
run = public.section_run()
scc = run.m_create(public.section_single_configuration_calculation)
system = run.m_create(public.section_system)
scc.single_configuration_calculation_to_system_ref = system

assert scc.single_configuration_calculation_to_system_ref == system

print(json.dumps(run.m_to_dict(), indent=2))
