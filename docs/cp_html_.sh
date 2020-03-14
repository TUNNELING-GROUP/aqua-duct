rm -rf $(find ../../TUNNELING-GROUP.github.io/aqua-duct/* | grep -v -F v0.)
cp -R build/html/* ../../TUNNELING-GROUP.github.io/aqua-duct/

CWD=`pwd`

cd ../../TUNNELING-GROUP.github.io/aqua-duct/

mv _static static
mv _images images
mv _modules modules
mv _sources sources

sed -i .qaz 's/_static/static/g' $(grep -Rl 'static' *)
sed -i .qaz 's/_images/images/g' $(grep -Rl 'images' *)
sed -i .qaz 's/_modules/modules/g' $(grep -Rl 'modules' *)
sed -i .qaz 's/_sources/sources/g' $(grep -Rl 'sources' *)

find ./ -name '*.qaz' -delete

git add .
git commit --message "docs update $(date)"
git push

cd $CWD

