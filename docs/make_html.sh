#!/bin/sh

# this builds current

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

#~ # installation pkg & RST file!
#~ rm -rf ../aquaduct*.tar.gz
#~ $( cd .. ; bash ./make_all_pkgs.sh )
#~ rm -rf source/aquaduct*.tar.gz
#~ mv ../aquaduct*.tar.gz source/
#~ cd source
#~ ls -1 aquaduct*.tar.gz | sort -r | awk '{print "* :download:`"$1"`"}' > aquaduct_download_list.rst
#~ cd ..

#~ # other versions docs
#~ echo 'Documentation for other versions of Aqua-Duct:\n' > source/other_versions.rst
#~ CURRENT='None'
#~ for tag in $( git tag | sort -r | head -1 )
#~ do
    #~ if [ $CURRENT = 'None' ]
    #~ then
        #~ CURRENT=$tag
        #~ echo '* `'$tag' <../current/index.html>`_ (current version)' >> source/other_versions.rst
    #~ else
        #~ echo '* `'$tag' <../'$tag'/index.html>`_' >> source/other_versions.rst
    #~ fi
#~ done
#~ echo '* `development version <../devel/index.html>`_ (use with care)' >> source/other_versions.rst

# AQ pip command
#AQPIP="install --extra-index-url https:\/\/testpypi.python.org\/pypi aqueduct"
AQPIP="install aquaduct"
sed -e 's/AQPIP/'"$AQPIP"'/' source/aquaduct_install.template > source/aquaduct_install.rst

# AQ installation requirements
echo "* Python 2.7 (CPython implementation)" > source/aquaduct_install_requires.rst
printf 'install_requires_nice(1)' | python2.7 -i ../src/setup.py -n --name | sed '1d' 2>&1 >> source/aquaduct_install_requires.rst | cat > /dev/null

rm -rf -- build/html*
$MAKE SPHINXBUILD=$SPHINXBUILD html

# rework links to other resources
#find build/html/ -iname '*.html' -exec sed -i 's/localhost/'$( hostname )'/g' {} +

rm -rf source/aquaduct*.tar.gz

export PYTHONPATH=$PYTHONPATH_CACHE

rm -rf aquaduct_docs.zip
( cd build/html ; zip -r -9 ../../aquaduct_docs.zip * )

#rsync -avz -P --delete build/html/ 192.168.1.15:/home/tljm/public_html/aq/
