#!/bin/bash

python -m celery -A nomad.processing worker -B -l info  -Q celery
