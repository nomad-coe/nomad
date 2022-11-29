#!/bin/bash

python -m nomad.cli admin ops gui-config
python -m uvicorn --host 0.0.0.0 nomad.app.main:app
