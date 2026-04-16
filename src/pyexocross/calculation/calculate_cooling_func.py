"""
Calculate cooling functions for molecular systems.

This module provides functions for calculating cooling functions which represent
the rate of energy loss through radiative transitions.
"""
import numpy as np
import numexpr as ne
from ..base.constants import hc, c2, PI
from ..gpu.base_gpu import using_gpu
from ..gpu.calculate_cooling_func_gpu import gpu_cal_cooling_func

# Calculate Cooling Function
def cal_cooling_func(A, v, Ep, gp, T, Q):
    """
    Calculate the cooling function at one or more temperatures.

    cooling_func = sum(A * h * c * v * gp * exp(-c2 * Ep / T)) / (4 * pi * Q)

    Parameters
    ----------
    A : np.ndarray
        Einstein A coefficient array, shape (n_levels,)
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    Ep : np.ndarray
        Upper state energy array, shape (n_levels,)
    gp : np.ndarray
        Upper state degeneracy array, shape (n_levels,)
    T : float
        Temperature in Kelvin
    Q : float
        Partition function at temperature T

    Returns
    -------
    float
        Cooling function value
    """
    if using_gpu():
        gpu_result = gpu_cal_cooling_func(A, v, Ep, gp, T, Q)
        if gpu_result is not None:
            return gpu_result

    _sum = ne.evaluate('sum(A * hc * v * gp * exp(-c2 * Ep / T))')  
    cf = ne.evaluate('_sum / (4 * PI * Q)')
    return cf
