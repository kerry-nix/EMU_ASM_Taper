# -*- coding: utf-8 -*-
# Copyright (c) St. Anne's University Hospital in Brno. International Clinical
# Research Center, Biomedical Engineering. All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

# Std imports

# Third pary imports

# Local imports


def try_jit_decorate(jit_kwargs):
    try:
        from numba import jit
        return jit(**jit_kwargs)
    except ImportError:
        return lambda x: x


def try_njit_decorate(jit_args, jit_kwargs):
    try:
        from numba import njit
        return njit(jit_args, **jit_kwargs)
    except ImportError:
        return lambda x: x
