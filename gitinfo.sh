#!/bin/sh
echo log: str = \'$(git log -1 --oneline)\'\\nref: str = \'$(git describe --all)\'\\nversion: str = \'$(git describe)\' > nomad/gitinfo.py