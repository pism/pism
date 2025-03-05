#!/bin/bash

# Replace paths to bundled JS libraries with paths to libraries installed using the package system

set -e
set -u

manual=$1
sphinxdoc="/usr/share/javascript/sphinxdoc/1.0"

for file in $(find ${manual} -name "*.html");
do
  for lib in jquery underscore doctools searchtools language_data;
  do
    sed -i "s%script src=\".*_static/\\(${lib}.js\\)\"%script src=\"${sphinxdoc}/\\1\"%" $file
    rm -f ${manual}/_static/${lib}.js
  done
done
