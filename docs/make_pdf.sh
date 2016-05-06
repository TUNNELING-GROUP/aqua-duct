#!/bin/sh

./make_html.sh

PYTHONPATH_CACHE=$PYTHONPATH
export PYTHONPATH=../src

if [ -n "`which gmake`" ]
then
    MAKE=gmake
else
    MAKE=make
fi

$MAKE latexpdf

export PYTHONPATH=$PYTHONPATH_CACHE

