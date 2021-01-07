#!/bin/bash
python -m nomad.cli admin ops gui-config
params=()
[ -e gunicorn.conf ] && params+=(--config gunicorn.conf)
[ -e gunicorn.log.conf ] && params+=(--log-config gunicorn.log.conf)
python -m gunicorn.app.wsgiapp "${params[@]}" -b 0.0.0.0:8000 nomad.app:app