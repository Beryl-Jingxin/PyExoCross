"""GPU tensor operations shared by spectroscopy calculations."""

from __future__ import annotations

from typing import Optional

import numpy as np

from pyexocross.base.constants import Inv8Pic, PI, c, c2, hc, hcInv4Pi
from pyexocross.gpu.base_gpu import _backend_arrays, _to_numpy, free_gpu_memory


_EPS = float(np.finfo(np.float64).eps)


def _safe_positive(mod, arr, eps: float = _EPS):
    return mod.where(arr > eps, arr, mod.full_like(arr, eps))


def _safe_nonzero(mod, arr, eps: float = _EPS):
    sign = mod.where(arr < 0, -1.0, 1.0)
    return mod.where(mod.abs(arr) > eps, arr, sign * eps)


def _finite_or_zero(mod, arr):
    return mod.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def gpu_cal_emicoefs(T_list, Q_list, Ep, gp, A, v, abundance) -> Optional[np.ndarray]:
    packed = _backend_arrays(T_list, Q_list, Ep, gp, A, v)
    if packed is None:
        return None
    provider, mod, arrs = packed
    T_g, Q_g, Ep_g, gp_g, A_g, v_g = arrs

    T_g = T_g[:, None]
    Q_g = Q_g[:, None]
    Ep_g = Ep_g[None, :]
    gp_g = gp_g[None, :]
    A_g = A_g[None, :]
    v_g = v_g[None, :]

    T_safe = _safe_positive(mod, T_g)
    Q_safe = _safe_nonzero(mod, Q_g)

    emiconst = gp_g * A_g * v_g * hcInv4Pi * abundance
    out = mod.exp(-c2 * Ep_g / T_safe) / Q_safe * emiconst
    out = _finite_or_zero(mod, out)
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result


def gpu_cal_emicoefs_nlte_2T(Tvib_list, Trot_list, Qnlte_arr, Evibp, Erotp, gp, A, v, abundance) -> Optional[np.ndarray]:
    packed = _backend_arrays(Tvib_list, Trot_list, Qnlte_arr, Evibp, Erotp, gp, A, v)
    if packed is None:
        return None
    provider, mod, arrs = packed
    Tvib_g, Trot_g, Q_g, Evibp_g, Erotp_g, gp_g, A_g, v_g = arrs

    Tvib_g = Tvib_g[:, None]
    Trot_g = Trot_g[:, None]
    Q_g = Q_g[:, None]
    Evibp_g = Evibp_g[None, :]
    Erotp_g = Erotp_g[None, :]
    gp_g = gp_g[None, :]
    A_g = A_g[None, :]
    v_g = v_g[None, :]

    Tvib_safe = _safe_positive(mod, Tvib_g)
    Trot_safe = _safe_positive(mod, Trot_g)
    Q_safe = _safe_nonzero(mod, Q_g)

    emiconst = gp_g * A_g * v_g * hcInv4Pi * abundance
    exponent = -c2 * (Evibp_g / Tvib_safe + Erotp_g / Trot_safe)
    out = mod.exp(exponent) / Q_safe * emiconst
    out = _finite_or_zero(mod, out)
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result


def gpu_cal_emicoefs_nlte_nvib(Trot_list, Qnlte_arr, nvib, Erotp, gp, A, v, abundance) -> Optional[np.ndarray]:
    packed = _backend_arrays(Trot_list, Qnlte_arr, nvib, Erotp, gp, A, v)
    if packed is None:
        return None
    provider, mod, arrs = packed
    Trot_g, Q_g, nvib_g, Erotp_g, gp_g, A_g, v_g = arrs

    Trot_g = Trot_g[:, None]
    Q_g = Q_g[:, None]
    nvib_g = nvib_g[None, :]
    Erotp_g = Erotp_g[None, :]
    gp_g = gp_g[None, :]
    A_g = A_g[None, :]
    v_g = v_g[None, :]

    Trot_safe = _safe_positive(mod, Trot_g)
    Q_safe = _safe_nonzero(mod, Q_g)

    emiconst = gp_g * A_g * v_g * hcInv4Pi * nvib_g * abundance
    out = mod.exp(-c2 * (Erotp_g / Trot_safe)) / Q_safe * emiconst
    out = _finite_or_zero(mod, out)
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result


def gpu_cal_emicoefs_nlte_pop(pop, A, v, abundance) -> Optional[np.ndarray]:
    packed = _backend_arrays(pop, A, v)
    if packed is None:
        return None
    provider, mod, arrs = packed
    pop_g, A_g, v_g = arrs

    out = (A_g * v_g * hcInv4Pi * pop_g * abundance)[None, :]
    out = _finite_or_zero(mod, out)
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result
