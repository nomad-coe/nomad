import sys

module_prefix = 'nomad.app.v1.models.graph'

setattr(sys.modules[f'{module_prefix}.utils'], 'ref_prefix', '#/definitions')
model = getattr(sys.modules[f'{module_prefix}.graph_models'], sys.argv[1])
print(model.schema_json(indent=2))
