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

    sum_val = mod.sum(A_g * hc * v_g * gp_g * mod.exp(-c2 * Ep_g / T_g))
    cf = sum_val / (4 * PI * Q_g)
    result = float(_to_numpy(provider, mod, cf))
    free_gpu_memory()
    return result
