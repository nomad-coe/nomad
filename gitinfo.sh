#!/bin/sh
echo log, ref, version, commit = \"$(git log -1 --oneline)\", \"$(git describe --all)\", \"$(git describe)\", \"$(git rev-parse --verify HEAD)\" > nomad/gitinfo.py