#!/bin/bash

python -m nomad.cli admin run app --with-gui --gunicorn --host 0.0.0.0 $@
