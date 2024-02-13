#!/bin/sh

set -e

working_dir=$(pwd)
project_dir=$(dirname $(dirname $(realpath $0)))

cd $project_dir/examples/data

rm -rf uploads/*.zip

zip -r -j uploads/theory.zip theory/*
zip -r -j uploads/eln.zip eln/*
zip -r -j uploads/tabular.zip tabular/*

cd cow_tutorial
zip -r ../uploads/cow_tutorial.zip *
cd ..

zip -r -j uploads/apm.zip apm/*
zip -r -j uploads/mpes.zip mpes/*
zip -r -j uploads/ellips.zip ellips/*
zip -r -j uploads/em_spctrscpy.zip em_spctrscpy/*
zip -r -j uploads/em_om.zip em_om/*
zip -r -j uploads/em_nion.zip em_nion/*
zip -r -j uploads/iv_temp.zip iv_temp/*
zip -r -j uploads/xps.zip xps/*
zip -r -j uploads/sts.zip sts/*
zip -r -j uploads/stm.zip stm/*
