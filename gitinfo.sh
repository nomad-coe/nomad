#!/bin/sh
echo log, ref, version, commit = \"$(git log -1 --oneline)\", \"$(git describe --all)\", \"$(git describe --tags)\", \"$(git rev-parse --verify --short HEAD)\" > nomad/gitinfo.py