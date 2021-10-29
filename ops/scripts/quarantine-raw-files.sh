#!/bin/sh
dir="quarantined"
pwd=`pwd`

[ -e raw-restricted.plain.zip ] && unzip -d $dir -o raw-restricted.plain.zip $@ > /dev/null
[ -e raw-public.plain.zip ] && unzip -d $dir -o raw-public.plain.zip $@ > /dev/null

if [ -e $dir ]
then
    cd $dir
    zip -rq $pwd/raw-quarantined.plain.zip *
    cd $pwd
    pwd
fi

[ -e raw-restricted.plain.zip ] && zip -dq raw-restricted.plain.zip $@
[ -e raw-public.plain.zip ] && zip -dq raw-public.plain.zip $@

[ -e $dir ] && rm -rf $dir

echo "done"
