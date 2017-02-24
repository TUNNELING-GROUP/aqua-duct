#!/usr/bin/env bash

CWD=`pwd`

RDIR=`mktemp -d`
TDIR=`mktemp -d`
PDIR=$CWD/build/current
rm -rf -- $PDIR
mkdir -p $PDIR

function make_docs {
    cd docs
    ./make_pdf.sh
    cd ..
    git stash
}

# prepare repo
cd $RDIR
git clone $CWD/.. aqua-duct
cd aqua-duct
git config user.email "tljm@wp.pl"
git config user.name "TM"


# make for highest tag
for tag in $( git tag | sort -r | head -1 )
#for tag in $( git branch | tail -n 1 )
do
    echo $tag
    git checkout tags/$tag # -b branch_$tag
    #git checkout $tag # -b branch_$tag
    VERSION=current
    #echo $tag $VERSION
    make_docs
    mv docs/build/html $PDIR
    cp docs/build/latex/*.pdf $PDIR
done

rm -rf -- $TDIR
rm -rf -- $RDIR
