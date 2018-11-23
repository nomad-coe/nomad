#!/bin/sh

cd .dependencies/nomad-lab-base
sbt docker repoTool
sbt docker repoWebservice