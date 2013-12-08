#!/bin/bash
for scriptname in $(ls p10km*sh) ; do
  echo
  echo "running ${scriptname} ..."
  bash $scriptname                       # will wait for completion
done
