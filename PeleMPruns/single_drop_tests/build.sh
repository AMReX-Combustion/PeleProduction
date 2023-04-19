#!/bin/sh

export PELE_PHYSICS_HOME=../../Submodules/PelePhysics
git pull
git submodule init
git submodule update
mkdir ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/multi_dechep
cp chem_files/* ${PELE_PHYSICS_HOME}/Support/Mechanism/Models/multi_dechep/
make -j8 TPL
make -j8
