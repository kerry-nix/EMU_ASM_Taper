# A MAkefile for the FSL installer

include ${FSLCONFDIR}/default.mk

PROJNAME = installer

pyinstall:
	fslpython -m pip install . --prefix ${FSLDEVDIR} --no-deps --ignore-installed --no-cache-dir -vvv
