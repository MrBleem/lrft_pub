#!/bin/bash

git config --global user.email "1194189468@qq.com"
git config --global user.name "MrBleem"

# echo "# lrft" >> README.md
version="v1.0.3"

git init
git add .
git commit -m "lrft_pub"
git branch -M main
git tag -a ${version} -m "try the function of release"
git remote add origin git@github.com:Project-XWH/lrft_pub.git
git push origin ${version}
# git pull origin v1
# git push -u origin v1


