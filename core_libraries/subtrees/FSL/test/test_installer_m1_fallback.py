#!/usr/bin/env python

import os
import os.path  as op
import contextlib
import shutil
import json

import fsl.installer.fslinstaller as inst

try:
    from unittest import mock
except ImportError:
    import mock

import pytest

from . import (server,
               CaptureStdout,
               indir,
               mock_input,
               strip_ansi_escape_sequences)


# mock miniconda installer which creates
# a mock $FSLDIR/bin/conda command
mock_miniconda_sh = """
#!/usr/bin/env bash

#called like <script> -b -p <prefix>
prefix=$3

mkdir -p $prefix/bin/
mkdir -p $prefix/etc/
mkdir -p $prefix/pkgs/

prefix=$(cd $prefix && pwd)

# called like
#  - conda env update -n base -f <envfile>
#  - conda clean -y --all
echo "#!/usr/bin/env bash"  >> $3/bin/conda
echo 'if   [ "$1" = "clean" ]; then '      >> $3/bin/conda
echo "    touch $prefix/cleaned"           >> $3/bin/conda
echo 'elif [ "$1" = "env" ]; then '        >> $3/bin/conda
echo "    cp "'$6'" $prefix/"              >> $3/bin/conda
echo "fi"                                  >> $3/bin/conda
chmod a+x $prefix/bin/conda
""".strip()


mock_manifest = """
{{
    "installer" : {{
        "version" : "1.3.4",
        "url"     : "na",
        "sha256"  : "na"
    }},
    "miniconda" : {{
        "macos-64" : {{
            "url"    : "{url}/miniconda.sh",
            "sha256" : "{condasha256}"
        }},
        "macos-M1" : {{
            "url"    : "{url}/miniconda.sh",
            "sha256" : "{condasha256}"
        }}
    }},
    "versions" : {{
        "latest" : "6.2.0",
        "6.2.0"  : [
            {{
                "platform"      : "macos-64",
                "environment"   : "{url}/env-macos-64-6.2.0.yml",
                "sha256"        : "{env64620sha256}",
                "base_packages" : ["fsl-base"],
                "output"        : {{
                    "install"   : {{ "version" : "2", "value" : "100" }}
                }}
            }},
            {{
                "platform"      : "macos-M1",
                "environment"   : "{url}/env-macos-M1-6.2.0.yml",
                "sha256"        : "{envM1620sha256}",
                "base_packages" : ["fsl-base"],
                "output"        : {{
                    "install"   : {{ "version" : "2", "value" : "100" }}
                }}
            }}
        ],
        "6.1.0"  : [
            {{
                "platform"      : "macos-64",
                "environment"   : "{url}/env-macos-64-6.1.0.yml",
                "sha256"        : "{env64610sha256}",
                "base_packages" : ["fsl-base"],
                "output"        : {{
                    "install"   : "100"
                }}
            }}
        ]
    }}
}}
""".strip()
# fmtvars: url, condasha256, env64620sha256, envM1620sha256, env64610sha256


mock_env_yml_template = """
{platform}
{version}
dependencies:
 - fsl-base 1234.0
""".strip()


@contextlib.contextmanager
def installer_server(cwd=None):
    if cwd is None:
        cwd = '.'
    cwd = op.abspath(cwd)

    with indir(cwd), server(cwd) as srv:
        with open('miniconda.sh', 'wt') as f:
            f.write(mock_miniconda_sh)
        for plat, ver in (('macos-64', '6.2.0'),
                          ('macos-M1', '6.2.0'),
                          ('macos-64', '6.1.0')):
            with open('env-{}-{}.yml'.format(plat, ver), 'wt') as f:
                f.write(mock_env_yml_template.format(platform=plat,
                                                     version=ver))

        condasha256    = inst.sha256('miniconda.sh')
        env64620sha256 = inst.sha256('env-macos-64-6.2.0.yml')
        envM1620sha256 = inst.sha256('env-macos-M1-6.2.0.yml')
        env64610sha256 = inst.sha256('env-macos-64-6.1.0.yml')

        manifest = mock_manifest.format(
            url=srv.url,
            condasha256=condasha256,
            env64620sha256=env64620sha256,
            envM1620sha256=envM1620sha256,
            env64610sha256=env64610sha256)

        with open('manifest.json', 'wt') as f:
            f.write(manifest)

        yield srv


def check_install(fsldir, platform, version):
    etc    = op.join(fsldir, 'etc')
    fslver = op.join(etc, 'fslversion')
    envyml = op.join(etc, 'env-{}-{}.yml'.format(platform, version))
    with open(fslver, 'rt') as f:
        assert f.read().strip() == version
    assert op.exists(envyml)


@contextlib.contextmanager
def pkgutil(retcode):
    with inst.tempdir(change_into=False) as td:
        exe = op.join(td, 'pkgutil')
        with open(exe, 'wt') as f:
            f.write('#!/usr/bin/env bash\n')
            f.write('exit {}\n'.format(retcode))
        os.chmod(exe, 0o755)
        path = op.pathsep.join((td, os.environ['PATH']))
        with mock.patch.dict(os.environ, PATH=path):
            yield


def test_installer_M1_install():
    with inst.tempdir():
        with installer_server() as srv:
            with mock.patch('fsl.installer.fslinstaller.FSL_RELEASE_MANIFEST',
                            '{}/manifest.json'.format(srv.url)):

                # macos-M1 build available for 6.2.0
                with mock.patch('platform.system',  return_value='darwin'), \
                     mock.patch('platform.machine', return_value='arm64'), \
                     inst.tempdir() as cwd:
                    inst.main(['--homedir', cwd, '--dest', 'fsl'])
                    check_install('fsl', 'macos-M1', '6.2.0')

                # macos-M1 build not available for 6.1.0 -
                # should fallback to macos-64
                with mock.patch('platform.system',  return_value='darwin'), \
                     mock.patch('platform.machine', return_value='arm64'), \
                     inst.tempdir() as cwd:
                    inst.main(['--homedir', cwd, '--dest', 'fsl', '-V', '6.1.0'])
                    check_install('fsl', 'macos-64', '6.1.0')


# installig an intel build on M1
def test_installer_M1_install_rosetta():
    with inst.tempdir():
        with installer_server() as srv:
            with mock.patch('fsl.installer.fslinstaller.FSL_RELEASE_MANIFEST',
                            '{}/manifest.json'.format(srv.url)), \
                mock.patch('platform.system',  return_value='darwin'), \
                mock.patch('platform.machine', return_value='arm64'):

                # pkgutil reports that rosetta is
                # enabled - installation should proceed
                with inst.tempdir() as cwd, pkgutil(0):
                    inst.main(['--homedir', cwd, '--dest', 'fsl',
                               '-V', '6.1.0'])
                    check_install('fsl', 'macos-64', '6.1.0')

                # pkgutil reports that rosetta is *not*
                # enabled - installation should be
                # aborted
                with inst.tempdir() as cwd, pkgutil(1):
                    with pytest.raises(SystemExit) as e:
                        inst.main(['--homedir', cwd, '--dest', 'fsl',
                                   '-V', '6.1.0'])
                    assert e.value.code != 0
                    assert not op.exists('fsl')
