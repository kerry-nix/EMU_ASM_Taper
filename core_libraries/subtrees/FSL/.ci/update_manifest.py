#!/usr/bin/env python
#
# This script is called when a new version of the fslinstaller.py script is
# tagged. It opens a merge request on the GitLab fsl/conda/manifest project,
# to update the FSL release manifest file so that it contains the new installer
# version string.
#


import os
import re

from fsl_ci        import (USERNAME,
                           EMAIL,
                           indir,
                           tempdir,
                           sprun)
from fsl_ci.gitlab import (gen_repository_url,
                           open_merge_request,
                           gen_branch_name)


MANIFEST_PATH = 'fsl/conda/manifest'

COMMIT_MSG = 'MNT: Update fslinstaller version to latest [{}] ' \
             'in FSL installer manifest'

MERGE_REQUEST_MSG = """
This MR was automatically opened as a result of a new tag being added
to the fsl/installer> project. It updates the installer version in
fsl-release.yml to the latest installer version.
""".strip()


def update_manifest(version):

    with open('fsl-release.yml', 'rt') as f:
        lines = list(f.readlines())

    pat     = r'installer:.*'
    updated = False

    for i, line in enumerate(lines):
        if re.match(pat, line):
            lines[i] = f'installer: {version}\n'
            updated = True
            break

    if not updated:
        for line in lines:
            print(line.rstrip())
        raise RuntimeError('Could not find installer section '
                           'in fsl-release.yml')

    with open('fsl-release.yml', 'wt') as f:
        for line in lines:
            f.write(line)


def checkout_and_update_manifest(server, token, tag, base_branch):

    manifest_url  = gen_repository_url(MANIFEST_PATH, server, token)
    branch        = f'mnt/{base_branch}/installer-{tag}'
    branch        = gen_branch_name(branch, MANIFEST_PATH, server, token)
    msg           = COMMIT_MSG.format(tag)

    with tempdir():
        sprun(f'git clone {manifest_url} manifest')
        with indir('manifest'):
            sprun(f'git config user.name  {USERNAME}')
            sprun(f'git config user.email {EMAIL}')
            sprun(f'git checkout -b {branch} origin/{base_branch}')
            update_manifest(tag)
            sprun( 'git add *')
            sprun(f'git commit -m "{msg}"')
            sprun(f'git push origin {branch}')

    return branch


def main(server=None, token=None, tag=None):

    if server is None: server = os.environ['CI_SERVER_URL']
    if token  is None: token  = os.environ['FSL_CI_API_TOKEN']
    if tag    is None: tag    = os.environ['CI_COMMIT_TAG']

    destination = 'master'
    branch = checkout_and_update_manifest(server, token, tag, destination)
    open_merge_request(MANIFEST_PATH,
                       branch,
                       MERGE_REQUEST_MSG,
                       server,
                       token,
                       destination)


if __name__ == '__main__':
    main()
