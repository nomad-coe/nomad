#!/bin/bash
while ! nc -z rabbitmq 5672; do sleep 3; done
while ! nc -z mongo 27017; do sleep 3; done
while ! nc -z elastic 9200; do sleep 3; done
python -m celery worker -l info -A nomad.processing -Q celery
