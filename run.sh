#!/bin/bash
find gui/build -type f | xargs -L1 bash -c 'sed "s_/fairdi/nomad/latest/gui_$1/gui_g" $2 > /tmp/temp_file; cp /tmp/temp_file $2;' -- $1
params=()
[ -e gunicorn.conf ] && params+=(--config gunicorn.conf)
[ -e gunicorn.log.conf ] && params+=(--log-config, gunicorn.log.conf)
python -m gunicorn.app.wsgiapp "${params[@]}" -b 0.0.0.0:8000 nomad.app:app