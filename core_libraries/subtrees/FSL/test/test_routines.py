#!/usr/bin/env python

import os
import os.path as op
import shutil
import subprocess as sp
import sys
import textwrap as tw

# py3
try:
    from unittest import mock
# py2
except ImportError:
    import mock

import pytest

from . import onpath, server

import fsl.installer.fslinstaller as inst


def test_identify_plaform():
    tests = [
        [('linux',  'x86_64'), 'linux-64'],
        [('darwin', 'x86_64'), 'macos-64'],
        [('darwin', 'arm64'),  'macos-M1'],
    ]

    for info, expected in tests:
        sys, cpu = info
        with mock.patch('platform.system', return_value=sys), \
             mock.patch('platform.machine', return_value=cpu):
            assert inst.identify_platform() == expected


def test_Version():
    assert inst.Version('1')       == inst.Version('1')
    assert inst.Version('1.2')     == inst.Version('1.2')
    assert inst.Version('1.2.3')   == inst.Version('1.2.3')
    assert inst.Version('1.2.3')   <  inst.Version('1.2.4')
    assert inst.Version('1.2.3')   >  inst.Version('1.2.2')
    assert inst.Version('1.2.3')   >  inst.Version('1.2')
    assert inst.Version('1.2.3')   <  inst.Version('1.2.3.0')
    assert inst.Version('1.2.3.0') >  inst.Version('1.2.3')


def test_get_admin_password():
    sudo = tw.dedent("""
    #!/usr/bin/env bash
    echo -n "Password: "
    read -e password
    if [ "$password" = "password" ]; then exit 0
    else exit 1
    fi
    """).strip()

    with inst.tempdir() as cwd:

        path = op.pathsep.join((cwd, os.environ['PATH']))

        with open('sudo', 'wt') as f:
            f.write(sudo)
        os.chmod('sudo', 0o755)

        # right password first time
        with mock.patch.dict(os.environ, PATH=path), \
             mock.patch('getpass.getpass', return_value='password'):
            assert inst.get_admin_password() == 'password'

        # wrong, then right
        returnvals = ['wrong', 'password']
        def getpass(*a):
            return returnvals.pop(0)
        with mock.patch.dict(os.environ, PATH=path), \
             mock.patch('getpass.getpass', getpass):
            assert inst.get_admin_password() == 'password'

        # wrong wrong wrong
        returnvals = ['wrong', 'bad', 'no']
        def getpass(*a):
            return returnvals.pop(0)
        with mock.patch.dict(os.environ, PATH=path), \
             mock.patch('getpass.getpass', getpass):
            with pytest.raises(Exception):
                inst.get_admin_password()


def test_download_file():

    with inst.tempdir() as cwd:
        with open('file', 'wt') as f:
            f.write('hello\n')
        with server(cwd) as srv:

            url = '{}/file'.format(srv.url)

            inst.download_file(url, 'copy')

            with open('copy', 'rt') as f:
                 assert f.read() == 'hello\n'

        # download_file should also work
        # with a path to a local file
        os.remove('copy')
        inst.download_file('file', 'copy')

        with open('copy', 'rt') as f:
            assert f.read() == 'hello\n'


def test_patch_file():

    content = tw.dedent("""
    line1
    line2
    line3
    line4
    """).strip()

    with inst.tempdir():

        with open('file', 'wt') as f:
            f.write(content)
        inst.patch_file('file', 'line2', 1, 'newline2')
        expect = content.replace('line2', 'newline2')
        with open('file', 'rt') as f:
            assert f.read().strip() == expect

        with open('file', 'wt') as f:
            f.write(content)
        inst.patch_file('file', 'line2', 2, 'newline2')
        expect = content.replace('line2\nline3', 'newline2')
        with open('file', 'rt') as f:
            assert f.read().strip() == expect

        with open('file', 'wt') as f:
            f.write(content)
        inst.patch_file('file', 'noline', 1, 'newline')
        expect = content + '\n\nnewline'
        with open('file', 'rt') as f:
            assert f.read().strip() == expect


def test_configure_shell():

    template = tw.dedent("""
    # FSL Setup
    FSLDIR={}
    PATH=${{FSLDIR}}/share/fsl/bin:${{PATH}}
    export FSLDIR PATH
    . ${{FSLDIR}}/etc/fslconf/fsl.sh
    """).strip()

    with inst.tempdir() as homedir:

        # no profile file exists
        inst.configure_shell('bash', homedir, '/fsl')
        with open('.bash_profile', 'rt') as f:
            assert f.read().strip() == template.format('/fsl')

        # existing profile with FSL config already present
        inst.configure_shell('bash', homedir, '/fsl_new')
        with open('.bash_profile', 'rt') as f:
            assert f.read().strip() == template.format('/fsl_new')

        # existing profile without FSL config
        config = tw.dedent("""
        line1
        line2
        line3
        """).strip()
        with open('.bash_profile', 'wt') as f:
            f.write(config)
        inst.configure_shell('bash', homedir, '/fsl')

        expect = config + '\n\n' + template.format('/fsl')
        with open('.bash_profile', 'rt') as f:
            assert f.read().strip() == expect


def test_configure_matlab():

    template = tw.dedent("""
    % FSL Setup
    setenv( 'FSLDIR', '{}' );
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    fsldir = getenv('FSLDIR');
    fsldirmpath = sprintf('%s/etc/matlab',fsldir);
    path(path, fsldirmpath);
    clear fsldir fsldirmpath;
    """).strip()

    with inst.tempdir() as homedir:

        os.makedirs(op.join('Documents', 'MATLAB'))
        startupm = op.join('Documents', 'MATLAB', 'startup.m')

        # no startup.m exists
        inst.configure_matlab(homedir, '/fsl')
        with open(startupm, 'rt') as f:
            assert f.read().strip() == template.format('/fsl')

        # existing startup.m with FSL config already present
        inst.configure_matlab(homedir, '/fsl_new')
        with open(startupm, 'rt') as f:
            assert f.read().strip() == template.format('/fsl_new')

        # existing startup.m without FSL config
        config = tw.dedent("""
        line1
        line2
        line3
        """).strip()
        with open(startupm, 'wt') as f:
            f.write(config)
        inst.configure_matlab(homedir, '/fsl')

        expect = config + '\n\n' + template.format('/fsl')
        with open(startupm, 'rt') as f:
            assert f.read().strip() == expect


def test_list_available_versions():
    manifest = {'versions' : {'latest' : '6.1.0',
                              '6.1.0' : [{'platform'    : 'linux-64',
                                          'environment' : 'http://env.yml'},
                                         {'platform'    : 'macos-64',
                                          'environment' : 'http://env.yml'},
                                         {'platform'    : 'linux-64',
                                          'cuda'        : '10.2',
                                          'environment' : 'http://env.yml'}],
                              '6.2.0' : [{'platform'    : 'linux-64',
                                          'environment' : 'http://env.yml'},
                                         {'platform'    : 'macos-64',
                                          'environment' : 'http://env.yml'},
                                         {'platform'    : 'linux-64',
                                          'cuda'        : '10.2',
                                          'environment' : 'http://env.yml'}]}}
    inst.list_available_versions(manifest)


def test_download_install_miniconda():
    class MockObject(object):
        pass

    installer = tw.dedent("""
    #!/usr/bin/env bash
    # is run like <miniconda.sh> -b -p <prefix>
    mkdir -p $3
    touch $3/installed
    """)

    def gen_manifest(platform, port, sha256):
        return {
            'miniconda' : {
                platform : {
                    'url'    : 'http://localhost:{}/remote.sh'.format(port),
                    'sha256' : sha256,
                }
            }
        }

    with inst.tempdir() as cwd:
        os.mkdir('remote')
        with open(op.join('remote', 'remote.sh'), 'wt') as f:
            f.write(installer)
        sha256 = inst.sha256(op.join('remote', 'remote.sh'))

        with server('remote') as srv:

            destdir                  = op.join(cwd, 'miniconda')
            ctx                      = MockObject()
            ctx.args                 = MockObject()
            ctx.destdir              = destdir
            ctx.need_admin           = False
            ctx.admin_password       = None
            ctx.args.no_checksum     = False
            ctx.args.skip_ssl_verify = False
            ctx.args.miniconda       = None
            ctx.platform             = 'linux'
            ctx.manifest             = gen_manifest('linux', srv.port, sha256)
            ctx.environment_channels = []
            ctx.run                  = lambda f, *a, **kwa: f(*a, **kwa)

            inst.download_miniconda(ctx)
            inst.install_miniconda(ctx)

            assert op.exists(destdir)
            assert op.exists(op.join(destdir, 'installed'))
            assert op.exists(op.join(destdir, '.condarc'))

            shutil.rmtree(destdir)

            ctx.manifest = gen_manifest('linux', srv.port, 'bad')
            with pytest.raises(Exception):
                inst.download_miniconda(ctx)

            ctx.args.no_checksum = True
            inst.download_miniconda(ctx)

            with open(op.join('remote', 'remote.sh'), 'wt') as f:
                f.write(installer)
            sha256 = inst.sha256(op.join('remote', 'remote.sh'))


def test_self_update():

    ver    = inst.__version__
    newver = '{}.0.0'.format(int(ver.split('.')[0]) + 1)
    oldver = '{}.0.0'.format(int(ver.split('.')[0]) - 1)

    new_installer = tw.dedent("""
    #!/usr/bin/env python
    from __future__ import print_function
    print("new version")
    """).strip()

    script_template = tw.dedent("""
    #!/usr/bin/env python
    from __future__ import print_function

    manifest = {{
        'installer' : {{
            'version' : '{ver}',
            'url'     : 'http://localhost:{port}/new_installer.py',
            'sha256'  : '{checksum}'
        }}
    }}

    import fslinstaller as inst
    inst.self_update(manifest, '.', {check_checksum})
    print('old version')
    """).strip()

    with inst.tempdir() as cwd:
        with server() as srv:

            shutil.copyfile(inst.__absfile__, 'fslinstaller.py')
            with open('new_installer.py', 'wt') as f:
                f.write(new_installer)
            checksum = inst.sha256('new_installer.py')

            # new version available
            with open('script.py', 'wt') as f:
                f.write(script_template.format(ver=newver,
                                               port=srv.port,
                                               checksum=checksum,
                                               check_checksum=True))
            got = sp.check_output([sys.executable, 'script.py'])
            assert got.decode('utf-8').strip() == 'new version'

            # new version available, bad checksum
            with open('script.py', 'wt') as f:
                f.write(script_template.format(ver=newver,
                                               port=srv.port,
                                               checksum='bad',
                                               check_checksum=True))
            got = sp.check_output([sys.executable, 'script.py'])
            # drop warning messages
            got = got.decode('utf-8').strip().split('\n')[-1]
            assert got == 'old version'

            # new version available, bad checksum, but skip checksum
            with open('script.py', 'wt') as f:
                f.write(script_template.format(ver=newver,
                                               port=srv.port,
                                               checksum='bad',
                                               check_checksum=False))
            got = sp.check_output([sys.executable, 'script.py'])
            got = got.decode('utf-8').strip().split('\n')[-1]
            assert got == 'new version'

            # same version available
            with open('script.py', 'wt') as f:
                f.write(script_template.format(ver=ver,
                                               port=srv.port,
                                               checksum=checksum,
                                               check_checksum=True))
            got = sp.check_output([sys.executable, 'script.py'])
            assert got.decode('utf-8').strip() == 'old version'

            # old version available
            with open('script.py', 'wt') as f:
                f.write(script_template.format(ver=oldver,
                                               port=srv.port,
                                               checksum=checksum,
                                               check_checksum=True))
            got = sp.check_output([sys.executable, 'script.py'])
            assert got.decode('utf-8').strip() == 'old version'
