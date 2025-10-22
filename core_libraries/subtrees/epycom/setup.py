# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

from setuptools import setup, find_packages

setup(name='epycom',
      version='0.0b4',
      install_requires=['numpy', 'scipy', 'pandas', 'scikit-learn', 'numba<0.54'],
      description='Package for EEG data processing and analysis',
      url='',
      author='FNUSA-ICRC, BME',
      author_email='jan.cimbalnik@fnusa.cz, jan.cimbalnik@mayo.edu',
      license='BSD 3.0',
      packages=find_packages(),
      keywords='EEG epilepsy signal processing',
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Intended Audience :: Healthcare Industry',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 3',
          'Operating System :: OS Independent',
          'Topic :: Scientific/Engineering :: Medical Science Apps.'],
      setup_requires=['pytest_runner'],
      tests_requires=['pytest', 'pytest-benchmark'],
      zip_safe=False)
