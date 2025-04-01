#!/usr/bin/env bash

./make_html.sh

PYTHON=python3

if [ -n "`which gmake`" ]
then
    MAKE=gmake
else
    MAKE=make
fi

export PYTHONPATH=$(pwd)/../src:$PYTHONPATH
echo "THIS IS THE PYTHONPATH:"
echo $PYTHONPATH

$MAKE latexpdf

cp build/latex/Aqua-Duct.pdf build/html

#rsync -avz -P --delete build/html/ 192.168.1.15:/home/tljm/public_html/aq/

