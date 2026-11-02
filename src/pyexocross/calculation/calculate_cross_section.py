"""
Calculate cross sections using various line profiles.

This module provides functions for calculating cross sections using
Doppler, Lorentzian, Voigt, and other line profiles.
"""
import numpy as np
import numexpr as ne
from scipy.special import roots_hermite
from ..base import constants as bin_constants
from ..gpu.base_gpu import using_gpu
from ..gpu.calculate_line_profile_gpu import gpu_cross_section
from ..calculation.calcualte_line_profile import (
    Doppler_profile,
    Lorentzian_profile,
    SciPyVoigt_profile,
    SciPyWofzVoigt_profile,
    HumlicekVoigt_profile,
    PseudoVoigt_profile,
    BinnedGaussian_profile,
    BinnedLorentzian_profile,
    BinnedVoigt_profile,
    BinnedVoigt_bnormq,
    BinnedVoigt_lorenz,
)

def _cross_section_gpu(kind, wn_grid, v, coef, cutoff, alpha=None, gamma=None, sigma=None, hV=None, eta=None):
    """Dispatch cross-section calculation to GPU backend implementation."""
    if not using_gpu():
        return None
    return gpu_cross_section(
        kind,
        wn_grid,
        v,
        coef,
        cutoff,
        alpha=alpha,
        gamma=gamma,
        sigma=sigma,
        hV=hV,
        eta=eta,
        bin_size2=bin_constants.binSize2,
        inv_bin_size_pi_half=bin_constants.InvbinSizePIhalf,
    )


# Calculate Cross Sections
def cross_section_Doppler(wn_grid, v, alpha, coef, cutoff):
    """
    Calculate cross sections using Doppler (Gaussian) line profile.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array, shape (n_points,)
    v : np.ndarray
        Transition wavenumber centers, shape (n_lines,)
    alpha : np.ndarray
        Doppler HWHM array, shape (n_lines,)
    coef : np.ndarray
        Line strength coefficients, shape (n_lines,)
    cutoff : float or str
        Wing cutoff value in cm^-1 or 'None' for no cutoff

    Returns
    -------
    np.ndarray
        Cross-section array, shape (n_points,)
    """  
    if using_gpu():
        gpu_xsec = _cross_section_gpu('doppler', wn_grid, v, coef, cutoff, alpha=alpha)
        if gpu_xsec is not None:
            return gpu_xsec

    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = 0
        end = len(wn_grid)
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            Doppler = Doppler_profile(dv, alpha)
            _xsec[idx] = coef @ Doppler 
    else:
        start = max(0,wn_grid.searchsorted(min(v)-cutoff)-1)
        end = min(wn_grid.searchsorted(max(v)+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            if filter.sum() > 0:
                _dv = dv[filter]
                _alpha = alpha[filter]
                _coef = coef[filter]
                Doppler = Doppler_profile(_dv, _alpha)
                _xsec[idx] = _coef @ Doppler  
    xsec[start:end] += _xsec
    return xsec

def cross_section_Lorentzian(wn_grid, v, gamma, coef, cutoff):
    """
    Calculate cross sections using Lorentzian line profile.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array, shape (n_points,)
    v : np.ndarray
        Transition wavenumber centers, shape (n_lines,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_lines,)
    coef : np.ndarray
        Line strength coefficients, shape (n_lines,)
    cutoff : float or str
        Wing cutoff value in cm^-1 or 'None' for no cutoff

    Returns
    -------
    np.ndarray
        Cross-section array, shape (n_points,)
    """  
    if using_gpu():
        gpu_xsec = _cross_section_gpu('lorentzian', wn_grid, v, coef, cutoff, gamma=gamma)
        if gpu_xsec is not None:
            return gpu_xsec

    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = 0
        end = len(wn_grid)
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            Lorentzian = Lorentzian_profile(dv, gamma)
            _xsec[idx] = coef @ Lorentzian
    else:
        start = max(0,wn_grid.searchsorted(min(v)-cutoff)-1)
        end = min(wn_grid.searchsorted(max(v)+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            if filter.sum() > 0:
                _dv = dv[filter]
                _gamma = gamma[filter]
                _coef = coef[filter]
                Lorentzian = Lorentzian_profile(_dv, _gamma)
                _xsec[idx] = _coef @ Lorentzian
    xsec[start:end] += _xsec   
    return xsec

def cross_section_SciPyVoigt(wn_grid, v, sigma, gamma, coef, cutoff):
    """
    Calculate cross sections using SciPy Voigt line profile.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array, shape (n_points,)
    v : np.ndarray
        Transition wavenumber centers, shape (n_lines,)
    sigma : np.ndarray
        Gaussian standard deviation array, shape (n_lines,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_lines,)
    coef : np.ndarray
        Line strength coefficients, shape (n_lines,)
    cutoff : float or str
        Wing cutoff value in cm^-1 or 'None' for no cutoff

    Returns
    -------
    np.ndarray
        Cross-section array, shape (n_points,)
    """  
    if using_gpu():
        gpu_xsec = _cross_section_gpu('scipy_voigt', wn_grid, v, coef, cutoff, sigma=sigma, gamma=gamma)
        if gpu_xsec is not None:
            return gpu_xsec
        
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = 0
        end = len(wn_grid)
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            SciPyVoigt = SciPyVoigt_profile(dv, sigma, gamma)
            _xsec[idx] = coef @ SciPyVoigt   
    else:
        start = max(0,wn_grid.searchsorted(min(v)-cutoff)-1)
        end = min(wn_grid.searchsorted(max(v)+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            if filter.sum() > 0:
                _dv = dv[filter]
                _sigma = sigma[filter]
                _gamma = gamma[filter]
                _coef = coef[filter]
                SciPyVoigt = SciPyVoigt_profile(_dv, _sigma, _gamma)
                _xsec[idx] = _coef @ SciPyVoigt   
    xsec[start:end] += _xsec
    return xsec

def cross_section_SciPyWofzVoigt(wn_grid, v, sigma, gamma, coef, cutoff):
    """
    Calculate cross sections using SciPy wofz Voigt line profile.

    Uses Faddeeva function (wofz) for Voigt profile calculation.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array, shape (n_points,)
    v : np.ndarray
        Transition wavenumber centers, shape (n_lines,)
    sigma : np.ndarray
        Gaussian standard deviation array, shape (n_lines,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_lines,)
    coef : np.ndarray
        Line strength coefficients, shape (n_lines,)
    cutoff : float or str
        Wing cutoff value in cm^-1 or 'None' for no cutoff

    Returns
    -------
    np.ndarray
        Cross-section array, shape (n_points,)
    """
    if using_gpu():
        gpu_xsec = _cross_section_gpu('scipy_wofz_voigt', wn_grid, v, coef, cutoff, sigma=sigma, gamma=gamma)
        if gpu_xsec is not None:
            return gpu_xsec
        
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = 0
        end = len(wn_grid)
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            SciPyWofzVoigt = SciPyWofzVoigt_profile(dv, sigma, gamma)
            _xsec[idx] = coef @ SciPyWofzVoigt
    else:
        start = max(0,wn_grid.searchsorted(min(v)-cutoff)-1)
        end = min(wn_grid.searchsorted(max(v)+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            if filter.sum() > 0:
                _dv = dv[filter]
                _sigma = sigma[filter]
                _gamma = gamma[filter]
                _coef = coef[filter]
                SciPyWofzVoigt = SciPyWofzVoigt_profile(_dv, _sigma, _gamma)
                _xsec[idx] = _coef @ SciPyWofzVoigt
    xsec[start:end] += _xsec
    return xsec

def cross_section_HumlicekVoigt(wn_grid, v, alpha, gamma, coef, cutoff):
    """
    Calculate cross sections using Humlicek Voigt line profile.

    Uses Humlicek's rational approximation for efficient Voigt profile calculation.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array, shape (n_points,)
    v : np.ndarray
        Transition wavenumber centers, shape (n_lines,)
    alpha : np.ndarray
        Doppler HWHM array, shape (n_lines,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_lines,)
    coef : np.ndarray
        Line strength coefficients, shape (n_lines,)
    cutoff : float or str
        Wing cutoff value in cm^-1 or 'None' for no cutoff

    Returns
    -------
    np.ndarray
        Cross-section array, shape (n_points,)
    """
    if using_gpu():
        gpu_xsec = _cross_section_gpu('humlicek_voigt', wn_grid, v, coef, cutoff, alpha=alpha, gamma=gamma)
        if gpu_xsec is not None:
            return gpu_xsec
        
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = 0
        end = len(wn_grid)
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            HumlicekVoigt = HumlicekVoigt_profile(dv, alpha, gamma)
            _xsec[idx] = coef @ HumlicekVoigt
    else:
        start = max(0,wn_grid.searchsorted(min(v)-cutoff)-1)
        end = min(wn_grid.searchsorted(max(v)+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            if filter.sum() > 0:
                _dv = dv[filter]
                _alpha = alpha[filter]
                _gamma = gamma[filter]
                _coef = coef[filter]
                HumlicekVoigt = HumlicekVoigt_profile(_dv, _alpha, _gamma)
                _xsec[idx] = _coef @ HumlicekVoigt
    xsec[start:end] += _xsec
    return xsec

def cross_section_PseudoVoigt(wn_grid, v, hV, eta, coef, cutoff):
    """
    Calculate cross sections using Pseudo-Voigt line profile.

    Pseudo-Voigt is a weighted sum of Gaussian and Lorentzian profiles.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array, shape (n_points,)
    v : np.ndarray
        Transition wavenumber centers, shape (n_lines,)
    hV : np.ndarray
        Effective Voigt HWHM shared by the Gaussian and Lorentzian
        components, shape (n_lines,)
    eta : np.ndarray
        Mixing parameter array (0 = pure Gaussian, 1 = pure Lorentzian), shape (n_lines,)
    coef : np.ndarray
        Line strength coefficients, shape (n_lines,)
    cutoff : float or str
        Wing cutoff value in cm^-1 or 'None' for no cutoff

    Returns
    -------
    np.ndarray
        Cross-section array, shape (n_points,)
    """
    if using_gpu():
        gpu_xsec = _cross_section_gpu(
            'pseudo_voigt',
            wn_grid,
            v,
            coef,
            cutoff,
            hV=hV,
            eta=eta,
        )
        if gpu_xsec is not None:
            return gpu_xsec

    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = 0
        end = len(wn_grid)
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            PseudoVoigt = PseudoVoigt_profile(dv, hV, eta)
            _xsec[idx] = coef @ PseudoVoigt
    else:
        start = max(0,wn_grid.searchsorted(min(v)-cutoff)-1)
        end = min(wn_grid.searchsorted(max(v)+cutoff),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            if filter.sum() > 0:
                _dv = dv[filter]
                _hV = hV[filter]
                _eta = eta[filter]
                _coef = coef[filter]
                PseudoVoigt = PseudoVoigt_profile(_dv, _hV, _eta)
                _xsec[idx] = _coef @ PseudoVoigt
    xsec[start:end] += _xsec   
    return xsec

def binnedIntegral(wn_grid, v, coef, cutoff, kind, alpha=None, gamma=None, sigma=None):
    """Average each line over bins, independently of output window and chunks."""
    from scipy.special import erf
    width = float(bin_constants.binSize2) / 2
    half = width / 2
    limit = None if cutoff in (None, 'None') else float(cutoff)
    result = np.zeros_like(wn_grid, dtype=float)
    if kind == 'binned_voigt':
        roots, weights = roots_hermite(20)
        shift = np.sqrt(2) * sigma[:, None] * roots[None, :]
        if limit is None:
            norm = np.full(len(v), np.pi * np.sqrt(np.pi))
        else:
            norm = np.sum(weights * (np.arctan((limit-shift)/gamma[:, None])
                                     - np.arctan((-limit-shift)/gamma[:, None])), axis=1)
    elif kind == 'binned_lorentzian':
        norm = np.pi if limit is None else 2 * np.arctan(limit/gamma)
    for i, centre in enumerate(wn_grid):
        dv = centre - v
        lower, upper = dv-half, dv+half
        if limit is not None:
            lower = np.maximum(lower, -limit)
            upper = np.minimum(upper, limit)
        valid = upper > lower
        if not np.any(valid):
            continue
        lo, hi = lower[valid], upper[valid]
        if kind == 'binned_gaussian':
            integral = (erf(np.sqrt(np.log(2))*hi/alpha[valid])
                        - erf(np.sqrt(np.log(2))*lo/alpha[valid])) / 2
        elif kind == 'binned_lorentzian':
            integral = (np.arctan(hi/gamma[valid])-np.arctan(lo/gamma[valid])) / (norm if np.isscalar(norm) else norm[valid])
        else:
            integral = np.sum(weights * (np.arctan((hi[:, None]-shift[valid])/gamma[valid, None])
                                        - np.arctan((lo[:, None]-shift[valid])/gamma[valid, None])), axis=1) / norm[valid]
        result[i] = coef[valid] @ integral / width
    return result


def cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff):
    gpu = _cross_section_gpu('binned_gaussian', wn_grid, v, coef, cutoff, alpha=alpha) if using_gpu() else None
    if gpu is not None:
        return gpu
    return binnedIntegral(wn_grid, v, coef, cutoff, 'binned_gaussian', alpha=alpha)


def cross_section_BinnedLorentzian(wn_grid, v, gamma, coef, cutoff):
    gpu = _cross_section_gpu('binned_lorentzian', wn_grid, v, coef, cutoff, gamma=gamma) if using_gpu() else None
    if gpu is not None:
        return gpu
    return binnedIntegral(wn_grid, v, coef, cutoff, 'binned_lorentzian', gamma=gamma)


def cross_section_BinnedVoigt(wn_grid, v, sigma, gamma, coef, cutoff):
    gpu = _cross_section_gpu('binned_voigt', wn_grid, v, coef, cutoff, sigma=sigma, gamma=gamma) if using_gpu() else None
    if gpu is not None:
        return gpu
    return binnedIntegral(wn_grid, v, coef, cutoff, 'binned_voigt', sigma=sigma, gamma=gamma)
