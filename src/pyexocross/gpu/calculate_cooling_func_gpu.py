"""GPU tensor operations shared by spectroscopy calculations."""

from __future__ import annotations

from typing import Optional

import numpy as np

from pyexocross.base.constants import Inv8Pic, PI, c, c2, hc, hcInv4Pi
from pyexocross.gpu.base_gpu import (
    _backend_arrays,
    _to_numpy,
    _torch_dtype,
    free_gpu_memory,
    get_torch_device,
)


_EPS = float(np.finfo(np.float64).eps)


def _safe_positive(mod, arr, eps: float = _EPS):
    return mod.where(arr > eps, arr, mod.full_like(arr, eps))


def _safe_nonzero(mod, arr, eps: float = _EPS):
    sign = mod.where(arr < 0, -1.0, 1.0)
    return mod.where(mod.abs(arr) > eps, arr, sign * eps)


def _finite_or_zero(mod, arr):
    return mod.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def gpu_cal_cooling_func(A, v, Ep, gp, T, Q) -> Optional[float]:
    """GPU cooling function scalar. Returns None if GPU not active."""
    packed = _backend_arrays(A, v, Ep, gp)
    if packed is None:
        return None

    provider, mod, arrs = packed
    A_g, v_g, Ep_g, gp_g = arrs

    if provider == 'cupy':
        T_g = mod.float64(T)
        Q_g = mod.float64(Q)
    else:
        device = get_torch_device()
        dtype = _torch_dtype(mod, device)
        T_g = mod.tensor(float(T), dtype=dtype, device=device)
        Q_g = mod.tensor(float(Q), dtype=dtype, device=device)

    T_safe = _safe_positive(mod, T_g)
    Q_safe = _safe_nonzero(mod, Q_g)
    term = A_g * hc * v_g * gp_g * mod.exp(-c2 * Ep_g / T_safe)
    term = _finite_or_zero(mod, term)
    sum_val = mod.sum(term)
    cf = _finite_or_zero(mod, sum_val / (4 * PI * Q_safe))
    result = float(_to_numpy(provider, mod, cf))
    free_gpu_memory()
    return result
