from __future__ import annotations

import json
from pathlib import Path


def _definitions_dir() -> Path:
    return Path(__file__).resolve().parents[2] / 'nomad' / 'layouts' / 'definitions'


def test_builtin_layout_definitions_have_required_fields():
    for definition_file in sorted(_definitions_dir().glob('*.json')):
        if definition_file.name == 'defaults.json':
            continue

        with definition_file.open() as f:
            payload = json.load(f)

        assert payload['id']
        assert payload['label']
        assert payload['overview']['type']
