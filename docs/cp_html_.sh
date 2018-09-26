rm -rf ../../tljm.github.io/aq/*
cp -R build/html/* ../../tljm.github.io/aq

CWD=`pwd`

cd ../../tljm.github.io/aq/

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

