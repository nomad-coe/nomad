#!/bin/bash

while ! nc -z app 8000; do sleep 3; done

python -m nomad.cli admin run hub
