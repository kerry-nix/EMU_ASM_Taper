#!/usr/bin/env bash

# Here we call the update_manifest.py script, which opens a merge request
# on the fsl/conda/manifest repository, to update the latest available
# installer version in the manifest.
#
# The update_manifest.py script uses functionality from fsl-ci-rules, so
# we need to install that before running the script.
python -m pip install --upgrade pip
python -m pip install git+https://git.fmrib.ox.ac.uk/fsl/fsl-ci-rules.git
python ./.ci/update_manifest.py
