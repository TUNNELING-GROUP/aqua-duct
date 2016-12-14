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
export PYTHONPATH=`pwd`/../src

$SPHINX_APIDOC -f -e -o source/ ../src/aqueduct/
sed -i '/undoc/d' source/*.*.rst

# ubuntu 14.04 helper
cp ../src/ubuntu_mdanalysis_install_helper.sh source/

# installation pkg & RST file!
rm -rf ../aqueduct*.tar.gz
$( cd .. ; ./make_pkg.sh )
rm -rf source/aqueduct*.tar.gz
mv ../aqueduct*.tar.gz source/
AQUEDUCT=$( basename $( ls source/aqueduct*.tar.gz ) )
sed 's/AQUEDUCT/'$AQUEDUCT'/' source/aqueduct_install.template > source/aqueduct_install.rst

rm -rf -- build/html*
$MAKE SPHINXBUILD=$SPHINXBUILD html

# rework links to other resources
#find build/html/ -iname '*.html' -exec sed -i 's/localhost/'$( hostname )'/g' {} +

rm -rf source/aqueduct*.tar.gz

export PYTHONPATH=$PYTHONPATH_CACHE

rm -rf aqueduct_docs.zip
$( cd build/html ; zip -r -9 ../../aqueduct_docs.zip * )

#rsync -avz -P --delete build/html/ 192.168.1.15:/home/tljm/public_html/aq/
