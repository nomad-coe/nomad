#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

set -x # echo on

mkdir tmp

python -m nomad.cli dev gui-artifacts --output-directory tmp

diff gui/src/searchQuantities.json tmp/searchQuantities.json
diff gui/src/metainfo.json tmp/metainfo.json
diff gui/src/parserMetadata.json tmp/parserMetadata.json
diff gui/src/toolkitMetadata.json tmp/toolkitMetadata.json
diff gui/src/unitsData.js tmp/unitsData.js
diff gui/src/exampleUploads.json tmp/exampleUploads.json
diff gui/src/northTools.json tmp/northTools.json

NOMAD_CONFIG=gui/tests/nomad.yaml python -m nomad.cli dev gui-config >tmp/env.js
diff gui/tests/env.js tmp/env.js

# cleanup
rm -rf tmp


