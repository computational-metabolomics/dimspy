#!/bin/sh

# http://prjemian.github.io/gh-pages.html
make clean
make html
cd _build/html
# touch .nojekyll
tar czf /tmp/html.tgz .
cd ../../../
# git checkout --orphan gh-pages
git checkout gh-pages
# git branch
# git rm -rf *
# git rm docs
# git rm -f .gitignore .travis.yml
tar xzf /tmp/html.tgz
git add .
git commit -a -m "Update html docs"
git push origin gh-pages
git checkout master
