import sys

setattr(sys.modules["nomad.app.v1.models.utils"], "ref_prefix", "#/definitions")
model = getattr(sys.modules["nomad.app.v1.graph_models"], sys.argv[1])
print(model.schema_json(indent=2))
