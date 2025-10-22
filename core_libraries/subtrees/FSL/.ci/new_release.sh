#!/usr/bin/env bash
#
# Deploy a new version of the fslinstaller script.  This script is run
# every time a new tag is added to the fsl/conda/installer gitlab
# repository.
#
# This simply involves copying the installer into the deployment directory,
# denoted by the $FSLINSTALLER_DEPLOY_DIRECTORY environment variable.

set -e


# Make sure that the installer version matches the tag
scriptver=$(cat fsl/installer/fslinstaller.py |
              grep "__version__ = "           |
              cut -d " " -f 3                 |
              tr -d "'")

if [ "$scriptver" != "$CI_COMMIT_TAG" ]; then
  echo "Version in fslinstaller.py does not match tag! $scriptver != $CI_COMMIT_TAG"
  exit 1
fi

cp fsl/installer/fslinstaller.py $FSLINSTALLER_DEPLOY_DIRECTORY/
