rm -rf ../../tunneling-group.github.io/aqua-duct/*
cp -R build/html/* ../../tunneling-group.github.io/aqua-duct

CWD=`pwd`

cd ../../tunneling-group.github.io/aqua-duct/

mv _static static
mv _images images
mv _modules modules
mv _sources sources

sed -i 's/_static/static/g' `grep -Rl '_static' *`
sed -i 's/_images/images/g' `grep -Rl '_images' *`
sed -i 's/_modules/modules/g' `grep -Rl '_modules' *`
sed -i 's/_sources/sources/g' `grep -Rl '_sources' *`

git add .
git commit --message "docs update $(date)"
git push

cd $CWD

