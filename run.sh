#!/bin/sh

find gui/build -type f -exec sed -i "s_/fairdi/nomad/latest/gui_$1/gui_g" {} \;
python -m gunicorn.app.wsgiapp --config gunicorn.conf --log-config gunicorn.log.conf -b 0.0.0.0:8000 nomad.app:app