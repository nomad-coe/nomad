from astroid import scoped_nodes
from astroid import MANAGER


def register(linter):
    # Needed for registering the plugin.
    pass


def transform(cls: scoped_nodes.ClassDef):
    if any(getattr(base, 'name', None) == 'Document' for base in cls.bases):
        cls.locals['objects'] = [scoped_nodes.FunctionDef('objects')]

MANAGER.register_transform(scoped_nodes.ClassDef, transform)
