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
export PYTHONPATH=~/.local/lib/python2.7/site-packages:`pwd`/../src

$SPHINX_APIDOC -f -e -o source/ ../src/aquaduct/
sed -i '/undoc/d' source/*.*.rst

# ubuntu 14.04 helper
cp ../src/ubuntu_mdanalysis_install_helper.sh source/

# installation pkg & RST file!
rm -rf ../aquaduct*.tar.gz
$( cd .. ; ./make_pkg.sh )
rm -rf source/aquaduct*.tar.gz
mv ../aquaduct*.tar.gz source/
AQUADUCT=$( basename $( ls source/aquaduct*.tar.gz ) )
sed 's/AQUADUCT/'$AQUADUCT'/' source/aquaduct_install.template > source/aquaduct_install.rst
AQPIP="pip install --extra-index-url https:\/\/testpypi.python.org\/pypi aqueduct"
#AQPIP="pip install aquaduct"
sed -i -e 's/AQPIP/'"$AQPIP"'/' source/aquaduct_install.rst
echo "* Python 2.7 (CPython implementation)" > source/aquaduct_install_requires.rst
printf 'install_requires_nice(1)' | python -i ../src/setup.py -n --name | sed '1d' 2>&1 >> source/aquaduct_install_requires.rst | cat > /dev/null

rm -rf -- build/html*
$MAKE SPHINXBUILD=$SPHINXBUILD html

# rework links to other resources
#find build/html/ -iname '*.html' -exec sed -i 's/localhost/'$( hostname )'/g' {} +

rm -rf source/aquaduct*.tar.gz

export PYTHONPATH=$PYTHONPATH_CACHE

rm -rf aquaduct_docs.zip
( cd build/html ; zip -r -9 ../../aquaduct_docs.zip * )

#rsync -avz -P --delete build/html/ 192.168.1.15:/home/tljm/public_html/aq/
