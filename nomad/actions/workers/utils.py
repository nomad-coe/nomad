from nomad.infrastructure import setup
from nomad.parsing.parsers import import_all_parsers


def worker_process_initializer():
    """
    Function that can be called to initialize the worker process during
    startup, so that the first tasks execute faster.
    """
    setup()
    import_all_parsers()
