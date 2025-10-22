#!/usr/bin/env python

import os.path as op

from setuptools import setup, find_namespace_packages


basedir = op.dirname(__file__)


version = {}
with open(op.join(basedir, 'fsl', 'installer', 'fslinstaller.py')) as f:
    for line in f:
        if line.startswith('__version__ = '):
            exec(line, version)
            break
version = version['__version__']


with open(op.join(basedir, 'README.md'), 'rt') as f:
    readme = f.read()


setup(
    name='fslinstaller',
    version=version,
    description='Scripts to install and update FSL',
    long_description=readme,
    long_description_content_type='text/markdown',
    url='https://git.fmrib.ox.ac.uk/fsl/conda/installer',
    author='Paul McCarthy',
    author_email='paul.mccarthy@ndcn.ox.ac.uk',
    license='Apache License Version 2.0',

    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],

    packages=['fsl', 'fsl.installer'],
    entry_points={
        'console_scripts' : [
            'fslinstaller.py = fsl.installer.fslinstaller:main',
        ]
    }
)
