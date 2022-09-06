#!/bin/sh

export PELE_PHYSICS_HOME=../../../Submodules/PelePhysics
git pull
git submodule init
git submodule update
make -j8 TPL
make -j8
