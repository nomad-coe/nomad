#!/bin/sh

echo "LOC with pygount (pip install pygount)"

echo "backend:       `pygount nomad/ -s py | awk '{print $1}' | paste -sd+ - | bc`"
echo "backend tests: `pygount tests/ -s py | awk '{print $1}' | paste -sd+ - | bc`"
echo "frontend:      `pygount gui/src -s js | awk '{print $1}' | paste -sd+ - | bc`"