#!/bin/bash
python -m nomad.cli admin ops gui-config
params=()
[ -e gunicorn.conf ] && params+=(--config gunicorn.conf)
[ -e gunicorn.log.conf ] && params+=(--log-config gunicorn.log.conf)
python -m gunicorn.app.wsgiapp "${params[@]}" --timeout 120 --keep-alive 5 --worker-class=uvicorn.workers.UvicornWorker -b 0.0.0.0:8000 nomad.app.main:app
