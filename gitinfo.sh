#!/bin/sh
echo log, ref, version = \'$(git log -1 --oneline)\', \'$(git describe --all)\', \'$(git describe)\' > nomad/gitinfo.py