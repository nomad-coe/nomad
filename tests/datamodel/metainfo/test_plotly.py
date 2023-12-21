import json
from tests.normalizing.conftest import run_processing


def test_plotly_snapshot(raw_files_function):
    directory = 'tests/data/datamodel/metainfo/plotly'
    mainfile = 'plotly.schema.archive.yaml'
    plotly_archive = run_processing(directory, mainfile)

    f = open('tests/data/datamodel/metainfo/plotly/snapshot.archive.json')
    snapshot = json.load(f)
    f.close()

    figures = plotly_archive.data['figures']
    snapshot_figures = snapshot['figures']
    for i in range(0, 4):
        assert json.dumps(figures[i].figure, sort_keys=True) == json.dumps(
            snapshot_figures[i], sort_keys=True
        )
