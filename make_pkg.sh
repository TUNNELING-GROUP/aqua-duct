#!/bin/sh

TDIR=`mktemp -d`
PDIR=$TDIR/aqueduct
mkdir -p $PDIR
cp -r src/* $PDIR
rm -rf $PDIR/testing

if [ `uname` = 'OpenBSD' ]
then
    find $PDIR -iname '*.pyc' -exec rm {} \;
    find $PDIR -iname '*.orig' -exec rm {} \;
else
    find $PDIR -iname '*.pyc' -delete
    find $PDIR -iname '*.orig' -delete
fi

VERSION=`echo '__version__' | python -i src/aqueduct/__init__.py 2>/dev/null | head -1 | cut -b 2- | rev | cut -b 2- | rev`
VERSION="_"$VERSION

tar -C $TDIR -cf aqueduct$VERSION.tar aqueduct
rm -rf $TDIR
gzip -f -9 aqueduct$VERSION.tar
