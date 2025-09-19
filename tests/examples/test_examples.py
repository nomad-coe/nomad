import os
import runpy

import pytest


@pytest.mark.parametrize(
    'path',
    [
        'examples/metainfo/data_frames.py',
        'examples/plugins',
    ],
)
def test_metainfo(path, capsys):
    """Runs the python files(s) in the given path."""
    abs_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../', path))
    if os.path.isdir(abs_path):
        files = find_py_files(abs_path)
    else:
        files = [abs_path]
    for file in files:
        runpy.run_path(file)

    capsys.readouterr()  # suppress stdout and stderr


def find_py_files(directory):
    """
    Recursively traverses the given directory and returns a list of absolute paths for all .py files.

    Args:
        directory (str): The path of the directory to traverse.

    Returns:
        list: A list of absolute paths to .py files.
    """
    # Iterate over all entries in the current directory
    py_files = []
    with os.scandir(directory) as entries:
        for entry in entries:
            # If it's a file and ends with '.py', add its absolute path
            if entry.is_file() and entry.name.endswith('.py'):
                py_files.append(os.path.abspath(entry.path))
            # If it's a directory, recursively search it
            elif entry.is_dir():
                py_files.extend(find_py_files(entry.path))

    return py_files
