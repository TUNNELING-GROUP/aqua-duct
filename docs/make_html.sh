#!/usr/bin/env bash

PYTHON=python2.7

# this builds current

SPHINX_APIDOC="sphinx-apidoc"
if [ -x ~/.local/bin/sphinx-apidoc ]
then
    SPHINX_APIDOC="$PYTHON $HOME/.local/bin/sphinx-apidoc"
fi
SPHINXBUILD="sphinx-build"
if [ -x ~/.local/bin/sphinx-build ]
then
    SPHINXBUILD="$PYTHON $HOME/.local/bin/sphinx-build"
fi
if [ -n "`which gmake`" ]
then
    MAKE=gmake
else
    MAKE=make
fi

#PYTHONPATH_CACHE=$PYTHONPATH
export PYTHONPATH=~/.local/lib/python2.7/site-packages:`pwd`/../src:$PYTHONPATH

$SPHINX_APIDOC -f -e -o source/ ../src/aquaduct/
sed -i '/undoc/d' source/*.*.rst

# ubuntu 14.04 helper
cp ../src/ubuntu_mdanalysis_install_helper.sh source/

# AQ pip command
#AQPIP="install --extra-index-url https:\/\/testpypi.python.org\/pypi aqueduct"
AQPIP="install aquaduct"
sed -e 's/AQPIP/'"$AQPIP"'/' source/aquaduct_install.template.rst > source/aquaduct_install.rst

# AQ installation requirements
echo "* Python 2.7 (CPython implementation)" > source/aquaduct_install_requires.rst
printf 'install_requires_nice(1)' | python2.7 -i ../src/setup.py -n --name 2> /dev/null > source/aquaduct_install_requires.rst

# valve HELP
sed '1,/HELP/!d' source/valve/valve_manual.rst.template | sed '$ d' > source/valve/valve_manual.rst
python2.7 ../src/apps/valve.py -h 2> /dev/null | awk '{print "    "$0}' >> source/valve/valve_manual.rst
sed '1,/HELP/d' source/valve/valve_manual.rst.template >> source/valve/valve_manual.rst

# pond HELP
sed '1,/HELP/!d' source/pond/pond_manual.rst.template | sed '$ d' > source/pond/pond_manual.rst
python2.7 ../src/apps/pond.py -h 2> /dev/null | awk '{print "    "$0}' >> source/pond/pond_manual.rst
sed '1,/HELP/d' source/pond/pond_manual.rst.template >> source/pond/pond_manual.rst

rm -rf -- build/html*
$MAKE SPHINXBUILD="$SPHINXBUILD" html 

# rework links to other resources
#find build/html/ -iname '*.html' -exec sed -i 's/localhost/'$( hostname )'/g' {} +

rm -rf source/aquaduct*.tar.gz

#export PYTHONPATH=$PYTHONPATH_CACHE

rm -rf aquaduct_docs.zip
#( cd build/html ; zip -r -9 ../../aquaduct_docs.zip * )

#rsync -avz -P --delete build/html/ 192.168.1.15:/home/tljm/public_html/aq/
