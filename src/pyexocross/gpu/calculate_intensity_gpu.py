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


def gpu_cal_abscoefs(T_list, Q_list, Epp, gp, A, v, abundance) -> Optional[np.ndarray]:
    packed = _backend_arrays(T_list, Q_list, Epp, gp, A, v)
    if packed is None:
        return None
    provider, mod, arrs = packed
    T_g, Q_g, Epp_g, gp_g, A_g, v_g = arrs

    T_g = T_g[:, None]
    Q_g = Q_g[:, None]
    Epp_g = Epp_g[None, :]
    gp_g = gp_g[None, :]
    A_g = A_g[None, :]
    v_g = v_g[None, :]

    T_safe = _safe_positive(mod, T_g)
    Q_safe = _safe_nonzero(mod, Q_g)
    v2_safe = _safe_nonzero(mod, v_g * v_g)

    absconst = gp_g * A_g * Inv8Pic / v2_safe * abundance
    out = mod.exp(-c2 * Epp_g / T_safe) * (1 - mod.exp(-c2 * v_g / T_safe)) / Q_safe * absconst
    out = _finite_or_zero(mod, out)
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result


def gpu_cal_abscoefs_nlte_2T(Tvib_list, Trot_list, Qnlte_arr, Evibpp, Erotpp, gp, A, v, abundance) -> Optional[np.ndarray]:
    packed = _backend_arrays(Tvib_list, Trot_list, Qnlte_arr, Evibpp, Erotpp, gp, A, v)
    if packed is None:
        return None
    provider, mod, arrs = packed
    Tvib_g, Trot_g, Q_g, Evibpp_g, Erotpp_g, gp_g, A_g, v_g = arrs

    Tvib_g = Tvib_g[:, None]
    Trot_g = Trot_g[:, None]
    Q_g = Q_g[:, None]
    Evibpp_g = Evibpp_g[None, :]
    Erotpp_g = Erotpp_g[None, :]
    gp_g = gp_g[None, :]
    A_g = A_g[None, :]
    v_g = v_g[None, :]

    Tvib_safe = _safe_positive(mod, Tvib_g)
    Trot_safe = _safe_positive(mod, Trot_g)
    Q_safe = _safe_nonzero(mod, Q_g)
    v2_safe = _safe_nonzero(mod, v_g * v_g)

    absconst = gp_g * A_g * Inv8Pic / v2_safe * abundance
    exponent = -c2 * (Evibpp_g / Tvib_safe + Erotpp_g / Trot_safe)
    out = mod.exp(exponent) * (1 - mod.exp(-c2 * v_g / Tvib_safe)) / Q_safe * absconst
    out = _finite_or_zero(mod, out)
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result


def gpu_cal_abscoefs_nlte_nvib(T_list, Trot_list, Qnlte_arr, nvib, Erotpp, gp, A, v, abundance) -> Optional[np.ndarray]:
    packed = _backend_arrays(T_list, Trot_list, Qnlte_arr, nvib, Erotpp, gp, A, v)
    if packed is None:
        return None
    provider, mod, arrs = packed
    T_g, Trot_g, Q_g, nvib_g, Erotpp_g, gp_g, A_g, v_g = arrs

    T_g = T_g[:, None]
    Trot_g = Trot_g[:, None]
    Q_g = Q_g[:, None]
    nvib_g = nvib_g[None, :]
    Erotpp_g = Erotpp_g[None, :]
    gp_g = gp_g[None, :]
    A_g = A_g[None, :]
    v_g = v_g[None, :]

    T_safe = _safe_positive(mod, T_g)
    Trot_safe = _safe_positive(mod, Trot_g)
    Q_safe = _safe_nonzero(mod, Q_g)
    v2_safe = _safe_nonzero(mod, v_g * v_g)

    absconst = gp_g * A_g * Inv8Pic / v2_safe * nvib_g * abundance
    exponent = -c2 * (Erotpp_g / Trot_safe)
    out = mod.exp(exponent) * (1 - mod.exp(-c2 * v_g / T_safe)) / Q_safe * absconst
    out = _finite_or_zero(mod, out)
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result


def gpu_cal_abscoefs_nlte_pop(T_list, pop, gp, gpp, A, v, abundance) -> Optional[np.ndarray]:
    packed = _backend_arrays(T_list, pop, gp, gpp, A, v)
    if packed is None:
        return None
    provider, mod, arrs = packed
    T_g, pop_g, gp_g, gpp_g, A_g, v_g = arrs

    T_g = T_g[:, None]
    pop_g = pop_g[None, :]
    gp_g = gp_g[None, :]
    gpp_g = gpp_g[None, :]
    A_g = A_g[None, :]
    v_g = v_g[None, :]

    T_safe = _safe_positive(mod, T_g)
    gpp_safe = _safe_nonzero(mod, gpp_g)
    v2_safe = _safe_nonzero(mod, v_g * v_g)

    absconst = gp_g / gpp_safe * A_g * Inv8Pic / v2_safe * pop_g * abundance
    out = (1 - mod.exp(-c2 * v_g / T_safe)) * absconst
    out = _finite_or_zero(mod, out)
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result
