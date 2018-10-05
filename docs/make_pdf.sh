#!/usr/bin/env bash

cd ~/Research/aqua-duct/docs


./make_html.sh

PYTHON=python2.7

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

export PYTHONPATH=~/.local/lib/python2.7/site-packages:`pwd`/../src:$PYTHONPATH

$MAKE SPHINXBUILD="$SPHINXBUILD" latexpdf

cp build/latex/Aqua-Duct.pdf build/html

rsync -avz -P --delete build/html/ 192.168.1.15:/home/tljm/public_html/aq/

