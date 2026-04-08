"""GPU tensor operations shared by spectroscopy calculations."""

from __future__ import annotations

from typing import Optional

import numpy as np

from pyexocross.base.constants import Inv8Pic, PI, c, c2, hc, hcInv4Pi
from pyexocross.gpu.base_gpu import _backend_arrays, _to_numpy, free_gpu_memory


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

    absconst = gp_g * A_g * Inv8Pic / (v_g ** 2) * abundance
    out = mod.exp(-c2 * Epp_g / T_g) * (1 - mod.exp(-c2 * v_g / T_g)) / Q_g * absconst
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

    absconst = gp_g * A_g * Inv8Pic / (v_g ** 2) * abundance
    exponent = -c2 * (Evibpp_g / Tvib_g + Erotpp_g / Trot_g)
    out = mod.exp(exponent) * (1 - mod.exp(-c2 * v_g / Tvib_g)) / Q_g * absconst
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

    absconst = gp_g * A_g * Inv8Pic / (v_g ** 2) * nvib_g * abundance
    exponent = -c2 * (Erotpp_g / Trot_g)
    out = mod.exp(exponent) * (1 - mod.exp(-c2 * v_g / T_g)) / Q_g * absconst
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

    absconst = gp_g / gpp_g * A_g * Inv8Pic / (v_g ** 2) * pop_g * abundance
    out = (1 - mod.exp(-c2 * v_g / T_g)) * absconst
    result = _to_numpy(provider, mod, out)
    free_gpu_memory()
    return result
