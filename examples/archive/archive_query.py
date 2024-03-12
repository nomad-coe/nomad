from nomad.client import ArchiveQuery
from nomad.metainfo import units

query = ArchiveQuery(
    query={
        'results.method.simulation.program_name': 'VASP',
        'results.material.elements': ['Ti', 'O'],
        'results.properties.geometry_optimization': {
            'final_energy_difference:lt': 1e-22,
        },
    },
    required={
        'run': '*',
        'workflow2': {
            'results': {
                'calculation_result_ref': {
                    'energy': '*',
                    'system_ref': {'chemical_composition_reduced': '*'},
                }
            }
        },
    },
)

for result in query.download(100):
    calc = None
    if result.workflow2.results:
        calc = result.workflow2.results.calculation_result_ref

    if not calc or not calc.system_ref:
        continue

    formula = calc.system_ref.chemical_composition_reduced
    if calc.energy.total:
        total_energy = calc.energy.total.value.to(units.eV)
    else:
        total_energy = 'N/A'
    print(f'{formula}: {total_energy}')

# Convert the retrieved data directly into a pandas dataframe
df_from_query = query.entries_to_dataframe(
    keys_to_filter=['workflow2.results.calculation_result_ref']
)
print(df_from_query)
