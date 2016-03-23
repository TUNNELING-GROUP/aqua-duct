#!/bin/sh

PYTHONPATH_CACHE=$PYTHONPATH
export PYTHONPATH=../src

sphinx-apidoc -f -e -o source/ ../src/aqueduct/
sed -i '/undoc/d' source/*.*.rst
gmake html

export PYTHONPATH=$PYTHONPATH_CACHE

