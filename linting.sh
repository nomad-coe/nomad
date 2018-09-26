#!/bin/sh

python -m pycodestyle --ignore=E501,E701 nomad tests
python -m pylint --load-plugins=pylint_mongoengine nomad tests
python -m mypy --ignore-missing-imports --follow-imports=silent --no-strict-optional nomad tests