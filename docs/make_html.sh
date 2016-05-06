#!/bin/sh

PYTHONPATH_CACHE=$PYTHONPATH
export PYTHONPATH=../src

sphinx-apidoc -f -e -o source/ ../src/aqueduct/
sed -i '/undoc/d' source/*.*.rst

if [ -n "`which gmake`" ]
then
    MAKE=gmake
else
    MAKE=make
fi

$MAKE html

export PYTHONPATH=$PYTHONPATH_CACHE

