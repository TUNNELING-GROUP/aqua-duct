#!/bin/sh

TDIR=`mktemp -d`
PDIR=$TDIR/aqueduct
mkdir -p $PDIR
cp -r src/* $PDIR
find $PDIR -iname '*.pyc' -delete
find $PDIR -iname '*.orig' -delete
tar -C $TDIR -cf aqueduct.tar aqueduct
rm -rf $TDIR
gzip -f -9 aqueduct.tar

