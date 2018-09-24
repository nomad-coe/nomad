#!/bin/sh
echo { \"log\": \"$(git log -1 --oneline)\", \"ref\": \"$(git describe --all)\", \"version\": \"$(git describe)\" } > src/gitinfo.json