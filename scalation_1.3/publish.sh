#!/usr/bin/env bash

CWD=$(pwd)

cd scalation_mathstat && sbt publishLocal && cd $CWD
cd scalation_modeling && sbt publishLocal && cd $CWD
cd scalation_models   && sbt publishLocal && cd $CWD

