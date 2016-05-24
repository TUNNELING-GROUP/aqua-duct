#!/bin/sh

SPHINX_APIDOC=sphinx-apidoc
if [ -x ~/.local/bin/sphinx-apidoc ]
then
    SPHINX_APIDOC=~/.local/bin/sphinx-apidoc
fi
SPHINXBUILD=sphinx-build
if [ -x ~/.local/bin/sphinx-build ]
then
    SPHINXBUILD=~/.local/bin/sphinx-build
fi
if [ -n "`which gmake`" ]
then
    MAKE=gmake
else
    MAKE=make
fi

PYTHONPATH_CACHE=$PYTHONPATH
export PYTHONPATH=../src

$SPHINX_APIDOC -f -e -o source/ ../src/aqueduct/
sed -i '/undoc/d' source/*.*.rst

# ubuntu 14.04 helper
cp ../src/ubuntu_mdanalysis_install_helper.sh source/
$( cd .. ; ./make_pkg.sh)
cp ../aqueduct.tar.gz source/

$MAKE SPHINXBUILD=$SPHINXBUILD html

# rework links to other resources
find build/html/ -iname '*.html' -exec sed -i 's/localhost/'$( hostname )'/g' {} +

export PYTHONPATH=$PYTHONPATH_CACHE

