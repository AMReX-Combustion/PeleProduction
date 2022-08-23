#!/bin/sh

export PELE_PHYSICS_HOME=../../../Submodules/PelePhysics
make -j8 TPL
make -j8
