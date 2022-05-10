#!/bin/sh
rm -rf uploads/*.zip
curl https://www.dropbox.com/s/8zd7aqe91lza2r4/theory-example-upload.zip?dl=0 -o uploads/theory.zip
# TODO this does not work on the mpcdf servers (no route to host) !?
# curl https://datashare.mpcdf.mpg.de/s/xeBEsWGyrRq9XH4/download -o uploads/theory.zip
cd eln
zip -r ../uploads/eln.zip *
cd ..
cd tabular
zip -r ../uploads/tabular.zip *
cd ..