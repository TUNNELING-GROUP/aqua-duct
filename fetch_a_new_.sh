

# target repo
REPO="test"
# init clean target reop
rm -rf $REPO $REPO.git
mkdir $REPO.git
cd $REPO.git
git init --bare
cd ..
git clone $REPO.git
cd $REPO
#git config --bool core.bare false
git config credential.helper store
cd ..

# source repo
SOURCE="source"
rm -rf $SOURCE
git clone https://gitlab.com/tljm/aqua-duct.git $SOURCE

cd $REPO

# GitLab updates ###############################################################

# import branches

for BRANCH in v0.2 v0.3 v0.4 v0.5 v0.5-no_exchange v1.0
do
	git fetch https://gitlab.com/tljm/aqua-duct.git $BRANCH
    if [ $? -eq 0 ]
    then
        git checkout FETCH_HEAD
        git checkout -b $BRANCH
        git push --set-upstream origin $BRANCH
    fi
done

# import tags
cd ../$SOURCE
TAGS=`git tag`
cd ../$REPO

for T in $TAGS
do
	git fetch https://gitlab.com/tljm/aqua-duct.git $T
	if [ $? -eq 0 ]
	then
		git checkout FETCH_HEAD
		git tag $T
		git push origin --tags
	fi
done

#exit 0

git remote set-url origin https://github.com/TUNNELING-GROUP/aqua-duct.git

# push to new remote
for BRANCH in v0.2 v0.3 v0.4 v0.5 v0.5-no_exchange v1.0
do
    git checkout $BRANCH
    git push
done
git push --tags

cd ..

