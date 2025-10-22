#!/usr/bin/env bash

set -e

git fetch origin master
thisver=$(  cat fsl/installer/fslinstaller.py |
            grep "__version__ = "             |
            cut -d " " -f 3                   |
            tr -d "'")
masterver=$(git show origin/master:fsl/installer/fslinstaller.py |
            grep "__version__ = "                                |
            cut -d " " -f 3                                      |
            tr -d "'")

if [ "$thisver" = "$masterver" ]; then
  echo "Version has not been updated!"
  echo "Version on master branch: $masterver"
  echo "Version in this MR:       $thisver"
  echo "The version number must be updated before this MR can be merged."
  exit 1
else
  echo "Version on master branch: $masterver"
  echo "Version in this MR:       $thisver"
fi
