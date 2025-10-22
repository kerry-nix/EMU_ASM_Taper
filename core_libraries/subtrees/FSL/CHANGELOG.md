# FSL installer script release history


# 3.5.0 (Wednesday 22nd March 2023)


 - Correctly determine the `root` user home directory in case the user
   has requested that the `root` user's shell profile should be modified.
 - New hidden `--debug` option, which enables very verbose output logging
   from `mamba` / `conda`.


# 3.4.2 (Sunday 12th March 2023)

 - Change the default installation directory to `/usr/local/fsl/` when the
   `fslinstaller.py` script is run as the root user. Additionally, do not
   modify the root user's shell profile.


# 3.4.1 (Wednesday 8th March 2023)

 - Make sure that the temporary installation directory is deleted as the root
   user if necessary.


# 3.4.0 (Thursday 2nd March 2023)

 - Fix the conda package cache directory (the `pkgs_dirs` setting) at
   `$FSLDIR/pkgs`, to avoid potential conflicts with user-configured package
   caches.
 - The installation log file is now copied to the user home directory on
   failure.


# 3.3.0 (Friday 27th January 2023)

 - Update the installer to install macOS-M1 FSL builds if available.
 - Exit with a warning if an Intel FSL build is to be installed on a
   M1 machine, and Rosetta emulation is not enabled.


# 3.2.1 (Tuesday 24th January 2023)

 - Unrecognised command-line arguments are ignored - this is to allow for
   forward-compatibility within a self-update cycle.
 - `bash` is used rather than `sh` when calling the miniconda installer
   script.


# 3.2.0 (Sunday 25th December 2022)

 - New hidden `--miniconda` option, allowing an alternate miniconda installer
   to be used.


# 3.1.0 (Saturday 24th December 2022)

 - Allow different progress reporting implementations
 - Clear all `$PYTHON*` environment variables before installing miniconda
   and FSL.


# 3.0.1 (Friday 13th December 2022)

 - Minor internal adjustments.


# 3.0.0

 - The installer script will now use `mamba` instead of `conda`, if present,
   for all conda commands.
 - Reverted to a single-step installation process - instead of installing
   base packages separately, the full installation is now performed with
   `conda env update -f <env>.yml`.
 - Use the number of package files saved to `$FSLDIR/pkgs/`to monitor
   and report progress of the main FSL installation, instead of counting
   the number of lines printed to standard output.

# 2.1.1

 - Added hooks to insert FSL license boilerplate into source files.

# 2.1.0

 - More internal changes and enhancements to improve usability in other
   scripts.


# 2.0.1

 - Internal changes to improve usability in other scripts.

# 2.0.0

 - Removed the `--cuda` / `--no_cuda` options.
 - Re-arrange the code to make it installable as a Python library.


# 1.10.2

 - Fix to handling of the `--cuda` / `--no_cuda` options on macOS.


# 1.10.1

 - Small adjustment to how the `devreleases.txt` file is parsed.

# 1.10.0

 - New hidden `--devrelease` and `--devlatest` options, for installing
   development releases.

# 1.9.0

 - Removed/disabled the `--update` option, for updating an existing FSL
   installation. This option may be re-enabled in the future.
 - Removed the hidden `--environment` option.
 - Update the `fslinstaller.py` script to work with the new CUDA package
   arrangement - FSL environment specifications are no longer provided
   for each supported CUDA version. Instead, all CUDA packages are included
   as part of the `linux-64` environment. The `--cuda` option can be used
   to select one set of packages to be installed, and the `--no_cuda` option
   can be used to exclude all CUDA packages from the installation.


# 1.8.0

 - The default FSL installation directory has been changed from `/usr/local/fsl/`
   to `$HOME/fsl`.
 - The fslinstaller now reads `FSLCONDA_USERNAME` and `FSLCONDA_PASSWORD` environment
   variables if a `--username` and `--password` were not supplied (only relevant for
   internal releases).
