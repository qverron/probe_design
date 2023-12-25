"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
import numpy as np  # type: ignore
import sys
from typing import Callable


def ert_type(x, stype, label):
    if not isinstance(x, stype):
        raise AssertionError(f"{label} should be {stype}, {type(x)} instead")


def ert_multiTypes(x, types, label):
    cond = any(isinstance(x, t) for t in types)
    if not cond:
        raise AssertionError(f"{label} should be one of {types}, {type(x)} instead")


def ert_nonNeg(x, label, include_zero=False):
    if not include_zero:
        if x <= 0:
            raise AssertionError(f"{label} should be greater than 0")
    elif x < 0:
        raise AssertionError(f"{label} should be greater than or equal to 0")


def ert_inInterv(x, vmin, vmax, label, leftClose=False, rightClose=True):
    if leftClose:
        if rightClose:
            if x < vmin or x > vmax:
                raise AssertionError(f"expected {vmin}<={label}<={vmax}")
        elif x < vmin or x >= vmax:
            raise AssertionError(f"expected {vmin}<={label}<{vmax}")
    elif rightClose:
        if x <= vmin or x > vmax:
            raise AssertionError(f"expected {vmin}<{label}<={vmax}")
    elif x <= vmin or x >= vmax:
        raise AssertionError(f"expected {vmin}<{label}<{vmax}")


def ert_in_dtype(x, dtype):
    if dtype.startswith("f"):
        if x > np.finfo(dtype).max:
            raise AssertionError(
                " ".join(
                    [
                        "expected to be lower than {dtype} max:",
                        f"{x} < {np.finfo(dtype).max}",
                    ]
                )
            )
    elif dtype.startswith("u") or dtype.startswith("i"):
        if x > np.iinfo(dtype).max:
            raise AssertionError(
                " ".join(
                    [
                        "expected to be lower than {dtype} max:",
                        f"{x} < {np.iinfo(dtype).max}",
                    ]
                )
            )
    else:
        logging.warning(f"assert not implemented for dtype '{dtype}'")


def enable_rich_assert(fun: Callable) -> Callable:
    def wrapper(*args, **kwargs):
        try:
            return fun(*args, **kwargs)
        except AssertionError as e:
            logging.exception(e)
            sys.exit()

    return wrapper
