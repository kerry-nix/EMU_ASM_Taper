#!/usr/bin/env python


import               os
import os.path    as op
import textwrap   as tw
import subprocess as sp

try:
    from unittest import mock
except ImportError:
    import mock

import pytest

from . import onpath, server

import fsl.installer.fslinstaller as inst


SUDO = """
#!/usr/bin/env bash
s=$1; shift
k=$2; shift

echo -n "Password: "
read password
echo $password > got_password

"$@"
""".strip()


def test_Process_check_call():
    with inst.tempdir() as cwd:
        script_template = tw.dedent("""
        #!/usr/bin/env sh
        touch {semaphore}
        exit {retcode}
        """).strip()

        with open('pass', 'wt') as f:
            f.write(script_template.format(semaphore='passed', retcode=0))
        with open('fail', 'wt') as f:
            f.write(script_template.format(semaphore='failed', retcode=1))
        os.chmod('pass', 0o755)
        os.chmod('fail', 0o755)

        inst.Process.check_call(op.join(cwd, 'pass'))
        assert op.exists('passed')

        with pytest.raises(Exception):
            inst.Process.check_call(op.join(cwd, 'fail'))
        assert op.exists('failed')


def test_Process_check_output():
    with inst.tempdir() as cwd:
        script_template = tw.dedent("""
        #!/usr/bin/env sh
        echo "{stdout}"
        exit {retcode}
        """).strip()

        # (stdout, retcode)
        tests = [
            ('stdout', 0),
            ('stdout', 1),
        ]

        for expect, retcode in tests:
            script = script_template.format(stdout=expect, retcode=retcode)

            with open('script', 'wt') as f:
                f.write(script)
            os.chmod('script', 0o755)

            if retcode == 0:
                got = inst.Process.check_output(op.join(cwd, 'script'))
                assert got.strip() == expect

            else:
                with pytest.raises(Exception):
                    inst.Process.check_output(op.join(cwd, 'script'))


def test_Process_monitor_progress():
    with inst.tempdir() as cwd:
        script = tw.dedent("""
        #!/usr/bin/env bash
        for ((i=0;i<10;i++)); do
            echo $i
        done
        touch $1
        """).strip()

        with open('script', 'wt') as f:
            f.write(script)

        os.chmod('script', 0o755)

        script  = op.join(cwd, 'script')

        # py2: make sure function accepts string and unicode
        scripts = [script, u'{}'.format(script)]

        for script in scripts:
            inst.Process.monitor_progress( script + ' a')
            inst.Process.monitor_progress([script + ' b'])
            inst.Process.monitor_progress([script + ' c', script + ' d'])
            inst.Process.monitor_progress( script + ' e',                 10)
            inst.Process.monitor_progress([script + ' f'],                10)
            inst.Process.monitor_progress([script + ' g', script + ' h'], 10)

            for touched in 'abcdefgh':
                assert op.exists(touched)
                os.remove(touched)

def test_Process_sudo_popen():
    with inst.tempdir() as cwd:
        cmd = tw.dedent("""
        #!/usr/bin/env bash

        echo "Running cmd" > command_output
        """)

        with open('sudo', 'wt') as f: f.write(SUDO)
        with open('cmd', 'wt')  as f: f.write(cmd)
        os.chmod('sudo', 0o755)
        os.chmod('cmd',  0o755)

        path = op.pathsep.join((cwd, os.environ['PATH']))

        with mock.patch.dict(os.environ, PATH=path):
            p = inst.Process.sudo_popen(['cmd'], 'password', stdin=sp.PIPE)
            p.communicate()

        with open('got_password', 'rt') as f:
            assert f.read().strip() ==  'password'
        with open('command_output', 'rt') as f:
            assert f.read().strip() == 'Running cmd'


def test_Process_popen_append_env():

    script = tw.dedent("""
    #!/usr/bin/env bash
    echo "$VAR1" >  output
    echo "$VAR2" >> output
    """).strip()

    with inst.tempdir():
        with open('script', 'wt') as f:
            f.write(script)
        os.chmod('script', 0o755)

        append = {'VAR1' : 'var1', 'VAR2' : 'var2'}

        p = inst.Process.popen('./script', append_env=append)
        p.wait()

        with open('output', 'rt') as f:
            assert f.read().strip() == 'var1\nvar2'


def test_Process_sudo_popen_append_env():

    script = tw.dedent("""
    #!/usr/bin/env bash
    echo "$VAR1" >  output
    echo "$VAR2" >> output
    """).strip()

    with inst.tempdir() as cwd:
        with open('sudo',   'wt') as f: f.write(SUDO)
        with open('script', 'wt') as f: f.write(script)
        os.chmod('script', 0o755)
        os.chmod('sudo',   0o755)

        append = {'VAR1' : 'var1', 'VAR2' : 'var2'}
        path   = op.pathsep.join((cwd, os.environ['PATH']))

        with mock.patch.dict(os.environ, PATH=path):
            p = inst.Process.sudo_popen(['script'],
                                        'password',
                                        stdin=sp.PIPE,
                                        append_env=append)
            p.wait()

        with open('output', 'rt') as f:
            assert f.read().strip() == 'var1\nvar2'
