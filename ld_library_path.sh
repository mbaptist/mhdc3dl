#!/bin/sh

PROJECTS_DIR=/home/mbaptist/fcup/projects/

CAT_LIBDIR=${PROJECTS_DIR}/cat/install/lib
GOOPS_LIBDIR=${PROJECTS_DIR}/goops/install/lib
LASS_LIBDIR=${PROJECTS_DIR}/lass/install/lib
GODESS_LIBDIR=${PROJECTS_DIR}/godess/install/lib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CAT_LIBDIR:$GOOPS_LIBDIR:$LASS_LIBDIR:$GODESS_LIBDIR

