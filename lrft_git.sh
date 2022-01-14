#!/bin/bash

git config --global user.email "1194189468@qq.com"
git config --global user.name "MrBleem"

# echo "# lrft" >> README.md
version="v1.1.1"

git init
git add .
git commit -m "lrft_pub"
git branch -M main
git tag -a ${version} -m "try the function of release"
git remote add origin git@github.com:Project-XWH/lrft_pub.git
git push origin ${version}
# # git pull origin v1
# # git push -u origin v1


# git filter-branch --index-filter 'git rm --cached --ignore-unmatch lrft.log'
# rm -rf .git/refs/original/
# git reflog expire --expire=now --all
# git fsck --full --unreachable
# git repack -A -d
# git gc --aggressive --prune=now
# git push --force origin ${version}