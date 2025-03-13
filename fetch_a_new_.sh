#!/bin/sh

# Git URLs
SOURCE_URL="git@gitlab.com:tljm/aqua-duct.git"
TARGET_URL="git@github.com:TUNNELING-GROUP/aqua-duct.git"

# Check if --push option is provided
PUSH=false
if [ "$1" = "--push" ]; then
    PUSH=true
else
    echo "Warning: --push option not used. Changes will not be pushed to the remote repository."
fi

# target repo
REPO="proto_target"
# init clean target repo
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
git clone $SOURCE_URL $SOURCE

cd $REPO

# GitLab updates ###############################################################

# import branches
cd ../$SOURCE
BRANCHES=$(git branch -r | grep -v '\->' | sed 's/ *origin\///')
echo "Branches to be imported: $BRANCHES"
cd ../$REPO

for BRANCH in $BRANCHES; do
    git fetch $SOURCE_URL $BRANCH
    if [ $? -eq 0 ]; then
        git checkout FETCH_HEAD
        git checkout -b $BRANCH
        git push --set-upstream origin $BRANCH
    fi
done

# import tags
cd ../$SOURCE
TAGS=$(git tag)
echo "Tags to be imported: $TAGS"
cd ../$REPO

for T in $TAGS; do
    git fetch $SOURCE_URL $T
    if [ $? -eq 0 ]; then
        git checkout FETCH_HEAD
        git tag $T
        git push origin --tags
    fi
done

# Ask user if they want to proceed with pushing to new remote
read -p "Do you want to proceed with pushing to new remote? (y/n): " proceed_remote
if [ "$proceed_remote" != "y" ]; then
    exit 0
fi

git remote set-url origin $TARGET_URL

# push to new remote
if [ "$PUSH" = true ]; then
    for BRANCH in $BRANCHES; do
        git checkout $BRANCH
        git push
    done
    git push --tags
else
    echo "Warning: Skipping push to new remote due to missing --push option."
fi

cd ..
