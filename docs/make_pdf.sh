#!/bin/sh

./make_html.sh

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

$MAKE SPHINXBUILD=$SPHINXBUILD latexpdf

export PYTHONPATH=$PYTHONPATH_CACHE

if [ -d ~/Dropbox/AQUADUCT/documentation_builds ]
then
    mkdir -p ~/Dropbox/AQUADUCT/documentation_builds/latex
    rm -rf ~/Dropbox/AQUADUCT/documentation_builds/latex/*
    cp -R build/latex/*.pdf ~/Dropbox/AQUADUCT/documentation_builds/latex
    cp -R build/html ~/Dropbox/AQUADUCT/documentation_builds
fi


