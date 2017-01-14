#!/bin/sh

CWD=`pwd`

RDIR=`mktemp -d`
TDIR=`mktemp -d`
PDIR=$TDIR/aquaduct
mkdir -p $PDIR

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

# make current
VERSION=current
make_docs
mv docs/build/html $PDIR/$VERSION

# make for each tag
for tag in $( git tag )
do
    git checkout tags/$tag
    VERSION=`make_version`
    make_docs
    mv docs/build/html $PDIR/$VERSION
done

#/home/tljm/dropbox_uploader.sh -p upload $CWD/*.tar.gz AQUADUCT/AQ_pkgs/

cd $TDIR
tar cvzf $CWD/AQ_full_help.tar.gz aquaduct

cd $CWD



rm -rf -- $TDIR
rm -rf -- $RDIR

