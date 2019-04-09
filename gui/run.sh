#!/bin/sh

find . -type f -exec sed -i "s_/fairdi/nomad/latest/gui_$1/gui_g" {} \;
nginx -g "daemon off;"