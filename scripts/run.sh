#!/bin/bash

python -m nomad.cli admin run app --with-gui --host 0.0.0.0 $@
