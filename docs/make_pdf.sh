#!/bin/sh

./make_html.sh

PYTHON=python2.7

SPHINX_APIDOC="$PYTHON sphinx-apidoc"
if [ -x ~/.local/bin/sphinx-apidoc ]
then
    SPHINX_APIDOC="$PYTHON $HOME/.local/bin/sphinx-apidoc"
fi
SPHINXBUILD="$PYTHON sphinx-build"
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

export PYTHONPATH=$PYTHONPATH:~/.local/lib/python2.7/site-packages:`pwd`/../src

$MAKE SPHINXBUILD="$SPHINXBUILD" latexpdf


