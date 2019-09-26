#!/bin/sh

set -e

git config -f .gitmodules --get-regexp '^submodule\..*\.path$' |
    while read path_key path
    do
        echo $path
        (echo "$path" | grep -vEq  ^dependencies/optimade-python-tools$) \
            && [ -f $path/requirements.txt ] && pip install -r $path/requirements.txt
        [ -f $path/setup.py ] && pip install --ignore-requires-python $1 $path
    done
