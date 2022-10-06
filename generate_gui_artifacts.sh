#!/bin/sh
python -m nomad.cli dev gui-artifacts --output-directory gui/src
python -m nomad.cli dev gui-config > gui/public/env.js
