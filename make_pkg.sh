#!/bin/sh

TDIR=`mktemp -d`
PDIR=$TDIR/aqueduct
mkdir -p $PDIR
cp -r src/* $PDIR

if [ `uname` = 'OpenBSD' ]
then
    find $PDIR -iname '*.pyc' -exec rm {} \;
    find $PDIR -iname '*.orig' -exec rm {} \;
else
    find $PDIR -iname '*.pyc' -delete
    find $PDIR -iname '*.orig' -delete
fi

tar -C $TDIR -cf aqueduct.tar aqueduct
rm -rf $TDIR
gzip -f -9 aqueduct.tar
