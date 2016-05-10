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

cp ../src/ubuntu_mdanalysis_install_helper.sh source/

$MAKE html

export PYTHONPATH=$PYTHONPATH_CACHE

