#!/usr/bin/env python
#
"""Utility functions used for testing. """


import os
import os.path as op
import contextlib
import threading
import multiprocessing as mp
import functools as ft
import sys
import re



# py3
try:
    import http.server as http
    from io import StringIO
    from unittest import mock

# py2
except ImportError:
    from StringIO import StringIO
    import SimpleHTTPServer as http
    http.HTTPServer = http.BaseHTTPServer.HTTPServer
    http.SimpleHTTPRequestHandler.protocol_version = 'HTTP/1.0'
    import mock


@contextlib.contextmanager
def onpath(dir):
    """Context manager which temporarily adds dir to the front of $PATH. """
    path               = os.environ['PATH']
    os.environ['PATH'] = dir + op.pathsep + path
    try:
        yield
    finally:
        os.environ['PATH'] = path


@contextlib.contextmanager
def indir(dir):
    """Context manager which temporarily changes into dir."""
    prevdir = os.getcwd()
    os.chdir(dir)
    try:
        yield
    finally:
        os.chdir(prevdir)


class HTTPServer(mp.Process):
    """Simple HTTP server which serves files from a specified directory.

    Intended to be used via the :func:`server` context manager function.
    """
    def __init__(self, rootdir):
        mp.Process.__init__(self)
        self.daemon = True
        self.rootdir = rootdir
        handler = http.SimpleHTTPRequestHandler
        self.server = http.HTTPServer(('', 0), handler)
        self.shutdown = mp.Event()

    def stop(self):
        self.shutdown.set()

    @property
    def port(self):
        return self.server.server_address[1]

    def run(self):
        with indir(self.rootdir):
            while not self.shutdown.is_set():
                self.server.handle_request()
            self.server.shutdown()


@contextlib.contextmanager
def server(rootdir=None):
    """Start a :class:`HTTPServer` on a separate thread to serve files from
    ``rootdir`` (defaults to the current working directory), then shut it down
    afterwards.
    """
    if rootdir is None:
        rootdir = os.getcwd()
    srv = HTTPServer(rootdir)
    srv.start()
    srv.url = 'http://localhost:{}'.format(srv.port)
    try:
        yield srv
    finally:
        srv.stop()


class CaptureStdout(object):
    """Context manager which captures stdout and stderr. """

    def __init__(self):
        self.reset()

    def reset(self):
        self.__mock_stdout = StringIO('')
        self.__mock_stderr = StringIO('')
        self.__mock_stdout.mode = 'w'
        self.__mock_stderr.mode = 'w'
        return self

    def __enter__(self):
        self.__real_stdout = sys.stdout
        self.__real_stderr = sys.stderr
        sys.stdout = self.__mock_stdout
        sys.stderr = self.__mock_stderr
        return self


    def __exit__(self, *args, **kwargs):
        sys.stdout = self.__real_stdout
        sys.stderr = self.__real_stderr
        return False

    @property
    def stdout(self):
        self.__mock_stdout.seek(0)
        return self.__mock_stdout.read()

    @property
    def stderr(self):
        self.__mock_stderr.seek(0)
        return self.__mock_stderr.read()


@contextlib.contextmanager
def mock_input(*responses):
    """Mock the built-in ``input`` or ``raw_input`` function so that it
    returns the specified sequence of ``responses``.

    Each response is returned from the ``input`` function in order, unless it
    is a callable, in which case it is called, and then the next non-callable
    response returned. This gives us a hacky way to manipulate things while
    stuck in an input REPL loop.
    """

    resp = iter(responses)

    def _input(*a, **kwa):
        n = next(resp)
        while callable(n):
            n()
            n = next(resp)
        return n

    if sys.version[0] == '2': target = '__builtin__.raw_input'
    else:                     target = 'builtins.input'

    with mock.patch(target, _input):
        yield


def strip_ansi_escape_sequences(text):
    """Does what function name says it does. """
    return re.sub(r'\x1b\[[0-9;]*m', '', text)
