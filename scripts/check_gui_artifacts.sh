#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir

set -x # echo on

mkdir tmp

python -m nomad.cli dev metainfo >tmp/metainfo.json
diff gui/src/metainfo.json tmp/metainfo.json

python -m nomad.cli dev search-quantities >tmp/searchQuantities.json
diff gui/src/searchQuantities.json tmp/searchQuantities.json

python -m nomad.cli dev toolkit-metadata >tmp/toolkitMetadata.json
diff gui/src/toolkitMetadata.json tmp/toolkitMetadata.json

python -m nomad.cli dev units >tmp/unitsData.js
diff gui/src/unitsData.js tmp/unitsData.js

python -m nomad.cli dev parser-metadata >tmp/parserMetadata.json
diff gui/src/parserMetadata.json tmp/parserMetadata.json

python -m nomad.cli dev gui-config >tmp/env.js
diff gui/public/env.js tmp/env.js

diff dependencies/nomad-remote-tools-hub/tools.json gui/src/northTools.json

python -m nomad.cli dev example-upload-metadata >tmp/exampleUploads.json
diff gui/src/exampleUploads.json tmp/exampleUploads.json

# cleanup
rm -rf tmp


