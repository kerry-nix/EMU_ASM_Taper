Core Libraries
===============

This folder is meant to house the scripts, wheels, and sutrees (etc.) that compose the core libraries used by the CNT.


# Folder Structure

## python

This folder stores the python wheels and environment files for the main research branches in the lab.

## matlab

Under development.

## Subtrees

Subtrees of external repositories. Used to manage a snapshot of our library and provide offline support, documentation, etc.

# Note

## Imaging

Many of the tools for imaging use dockerized containers or external libraries meant for ease of viewing. Where possible, we provide python and matlab support. As needed, we also provide external libraries and build files in other languages to aid in this research.

### FSL
A symbolic link to the FSL installer can be found at python/imaging/. As a note, FSL installs to a separate python environment. 

This environment can be activated using:

> conda activate fsal

## Dockers

### FSL
TBD

### FastSurfer
For more information, please visit: [FastSurfer Docker](https://github.com/Deep-MI/FastSurfer/tree/stable/Docker)
