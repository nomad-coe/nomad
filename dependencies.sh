#!/bin/sh

set -e

git config -f .gitmodules --get-regexp '^submodule\..*\.path$' |
    while read path_key path
    do
        echo $path
        [ -f $path/requirements.txt ] && pip install -r $path/requirements.txt
        [ -f $path/setup.py ] && pip install $1 $path
    done
