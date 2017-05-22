#!/usr/bin/env bash

mkdir -p data
echo "CLI options: $*"
sbt "run-main apps.GeneAnalysis2 $*"

