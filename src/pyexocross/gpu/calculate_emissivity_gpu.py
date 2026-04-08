"""GPU tensor operations shared by spectroscopy calculations."""

from __future__ import annotations

from typing import Optional

import numpy as np

from pyexocross.base.constants import Inv8Pic, PI, c, c2, hc, hcInv4Pi
from pyexocross.gpu.base_gpu import _backend_arrays, _to_numpy, free_gpu_memory


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

    emiconst = gp_g * A_g * v_g * hcInv4Pi * abundance
    out = mod.exp(-c2 * Ep_g / T_g) / Q_g * emiconst
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

    emiconst = gp_g * A_g * v_g * hcInv4Pi * abundance
    exponent = -c2 * (Evibp_g / Tvib_g + Erotp_g / Trot_g)
    out = mod.exp(exponent) / Q_g * emiconst
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

    emiconst = gp_g * A_g * v_g * hcInv4Pi * nvib_g * abundance
    out = mod.exp(-c2 * (Erotp_g / Trot_g)) / Q_g * emiconst
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
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result
