#!/bin/sh

CWD=`pwd`

gmake clean
./make_current.sh
./make_pdf.sh

TDIR=`mktemp -d`
PDIR=$TDIR/aquaduct
mkdir $PDIR

cp -R build/current/html/* $PDIR
mkdir $PDIR/devel
cp -R build/html/* $PDIR/devel

cd $PDIR
zip -r -9 $CWD/aquaduct_docs.zip *
cd $CWD
rm -rf -- $TDIR
rm -rf -- build/full
mkdir -p build/full
cd build/full
unzip $CWD/aquaduct_docs.zip
cd $CWD

