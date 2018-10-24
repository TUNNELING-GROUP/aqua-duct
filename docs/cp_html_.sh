#!/bin/sh

rm -rf ../../TUNNELING-GROUP.github.io/aqua-duct/*
cp -R build/html/* ../../TUNNELING-GROUP.github.io/aqua-duct

CWD=`pwd`

cd ../../TUNNELING-GROUP.github.io/aqua-duct/

mv _static static
mv _images images
mv _modules modules
mv _sources sources

sed -i 's/_static/static/g' `grep -Rl 'static' *`
sed -i 's/_images/images/g' `grep -Rl 'images' *`
sed -i 's/_modules/modules/g' `grep -Rl 'modules' *`
sed -i 's/_sources/sources/g' `grep -Rl 'sources' *`

git add .
git commit --message "docs update $(date)"
git push

cd $CWD

