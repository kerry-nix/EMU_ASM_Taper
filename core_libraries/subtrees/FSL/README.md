# FSL installer


[![pipeline status](https://git.fmrib.ox.ac.uk/fsl/conda/installer/badges/master/pipeline.svg)](https://git.fmrib.ox.ac.uk/fsl/conda/installer/-/commits/master)
[![coverage report](https://git.fmrib.ox.ac.uk/fsl/conda/installer/badges/master/coverage.svg)](https://git.fmrib.ox.ac.uk/fsl/conda/installer/-/commits/master)


This repository is the home of `fslinstaller.py`, the installer script for
[FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/).


The `fslinstaller.py` script in this repository is the successor to the
`fslinstaller.py` script from the fsl/installer> repository.  _This_ version
is for **conda-based** FSL release, from FSL version 6.0.6 onwards.


`fslinstaller.py`  is a Python script which can run with any version of
Python from 2.7 onwards. Normal usage of `fslinstaller.py` will look like
one of the following:

    ```
    python  fslinstaller.py
    python2 fslinstaller.py
    python3 fslinstaller.py
    curl https://some_url/fslinstaller.py | python
    ```


The `fslinstaller` script is also built as a Python package, and importable
via the `fsl.installer` package.  The conda package is called `fsl-installer`,
and is built at the fsl/conda/fsl-installer> repository.


In normal usage, the `fslinstaller.py` script performs the following tasks:
 1. Downloads the FSL release manifest file from a hard-coded URL, which is a
    JSON file containing information about available FSL releases.
 2. Asks the user where they would like to install FSL (hereafter referred to
    as `$FSLDIR`).
 3. Asks the user for their administrator password if necessary.
 4. Downloads a `miniconda` installer.
 5. Installs `miniconda` to `$FSLDIR`.
 6. Downloads a YAML file containing a conda environment specification for
    the latest FSL version (or the version requested by the user; hereafter
    referred to as `environment.yml`).
 7. Installs the FSL environment by running:
       `$FSLDIR/bin/conda env update -n base -f environment.yml`
 8. Modifies the user's shell configuration so that FSL is accessible in
    their shell environment.


Several advanced options are available - run `python fslinstaller.py -h`, and
read the `parse_args` function for more details on the advanced/hidden options.


# Managing `fslinstaller.py` versions and releases


All releases of `fslinstaller.py` are given a version of the form
`major.minor.patch`, for example `1.3.2`.

The fsl/conda/installer> project follows semantic versioning conventions,
where:
 - changes to the command-line interface require the major version number
   to be incremented
 - enhancements and new features require the minor version number to be
   incremented
 - bug fixes and minor changes require the patch version number to be
   incremented.

All changes to the `fslinstaller.py` must be accompanied by a change to the
`__version__` attribute in the `fslinstaller.py` script.


New versions of the `fslinstaller.py` script can be released simply by
creating a new tag, containing the new version identifier, on the
fsl/conda/installer> GitLab repository. This will cause the following
automated routines to run:

 - The new version of the `fslinstaller.py` script is deployed to a web server,
   available for download.

 - A merge request is opened on the fsl/conda/fsl-installer> conda recipe
   repository, causing the new version to be built as a conda package.

 - A merge request is optionally opened on the fsl/conda/manifest> repository,
   updating the installer version number in the FSL release manifest JSON
   file.


Note that the tag must be identical to the value of the `__version__`
attribute in the `fslinstaller.py` script.
