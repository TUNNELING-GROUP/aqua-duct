#!/bin/sh

CWD=`pwd`

RDIR=`mktemp -d`
TDIR=`mktemp -d`
PDIR=$TDIR/aquaduct
mkdir -p $PDIR

function copy2dir {
    cp -r src/README $PDIR
    cp -r src/apps $PDIR
    cp -r src/aquaduct $PDIR
    cp -r src/license.txt $PDIR
    cp -r src/setup.py $PDIR
    cp -r src/ubuntu_mdanalysis_install_helper.sh $PDIR

    rm -rf $PDIR/testing

    if [ `uname` = 'OpenBSD' ]
    then
        find $PDIR -iname '*.pyc' -exec rm {} \;
        find $PDIR -iname '*.orig' -exec rm {} \;
    else
        find $PDIR -iname '*.pyc' -delete
        find $PDIR -iname '*.orig' -delete
    fi
}

function make_version {
    VER=`echo '__version__' | python -i src/aquaduct/__init__.py 2>/dev/null | head -1 | cut -b 2- | rev | cut -b 2- | rev`
    echo "_"$VER
}

function make_pkg {
    rm -rf -- $CWD/aquaduct$1.tar
    tar -C $TDIR -cf $CWD/aquaduct$1.tar aquaduct
    rm -rf -- $CWD/aquaduct$1.tar.gz
    gzip -f -9 $CWD/aquaduct$1.tar
}

# prepare repo
cd $RDIR
git clone $CWD aqua-duct
cd aqua-duct

# make current
VERSION=_current
copy2dir
make_pkg $VERSION

# make for each tag
for tag in $( git tag )
do
    git checkout tags/$tag
    VERSION=`make_version`
    copy2dir
    make_pkg $VERSION
done

#/home/tljm/dropbox_uploader.sh -p upload $CWD/*.tar.gz AQUADUCT/AQ_pkgs/

cd $CWD

rm -rf -- $TDIR
rm -rf -- $RDIR
