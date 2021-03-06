#!/bin/bash

CWD=`pwd`

RDIR=`mktemp -d`
TDIR=`mktemp -d`
PDIR=$TDIR/aquaduct
#mkdir -p $PDIR

function make_version {
    VER=`echo '__version__' | python -i src/aquaduct/__init__.py 2>/dev/null | head -1 | cut -b 2- | rev | cut -b 2- | rev`
    echo $VER
}

function make_docs {
    cd docs
    pwd
    ./make_html.sh
    cd ..
    git stash
}

# prepare repo
cd $RDIR
git clone $CWD aqua-duct
cd aqua-duct
git config user.email "tljm@wp.pl"
git config user.name "TM"


# make for each tag
for tag in $( git tag | sort -r | head -1 )
do
    git checkout tags/$tag # -b branch_$tag
    VERSION=current
    #echo $tag $VERSION
    make_docs
    mv docs/build/html $PDIR
done

# make current
git checkout master
VERSION=devel
make_docs
mv docs/build/html $PDIR/$VERSION


#/home/tljm/dropbox_uploader.sh -p upload $CWD/*.tar.gz AQUADUCT/AQ_pkgs/

cd $TDIR
tar cvzf $CWD/AQ_full_help.tar.gz aquaduct
cd aquaduct
zip -r -9 $CWD/aquaduct_docs.zip *

cd $CWD
mkdir -p docs/build/full
cd docs/build/full
rm -rf -- *
tar xvzf ../../../AQ_full_help.tar.gz
cd $CWD


rm -rf -- $TDIR
rm -rf -- $RDIR
