#!/bin/sh

set -e
working_dir=`pwd`
echo $working_dir

git config -f .gitmodules --get-regexp '^submodule\..*\.path$' |
    while read path_key path
    do
        cd $working_dir
        cd $path
        if [ -z "$(git status --porcelain)" ]; then
            echo "$path is clean"
        else
            echo "$path is not clean"
            git stash
            git checkout nomad-fair-metainfo
            git stash apply
            git add -A
            git commit -a -m "Fixed gitignore."
            git push origin nomad-fair-metainfo
        fi
    done
