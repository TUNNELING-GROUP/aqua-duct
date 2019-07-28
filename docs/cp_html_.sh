rm -rf $(find ../../TUNNELING-GROUP.github.io/aqua-duct/ | grep -v -F v0.)
cp -R build/html/* ../../TUNNELING-GROUP.github.io/aqua-duct/

CWD=`pwd`

cd ../../TUNNELING-GROUP.github.io/aqua-duct/

mv static static
mv images images
mv modules modules
mv sources sources

sed -i 's/static/static/g' `grep -Rl 'static' *`
sed -i 's/images/images/g' `grep -Rl 'images' *`
sed -i 's/modules/modules/g' `grep -Rl 'modules' *`
sed -i 's/sources/sources/g' `grep -Rl 'sources' *`

git add .
git commit --message "docs update $(date)"
git push

cd $CWD

