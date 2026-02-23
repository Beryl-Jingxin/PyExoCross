"""
Calculate cross sections using various line profiles.

This module provides functions for calculating cross sections using
Doppler, Lorentzian, Voigt, and other line profiles.
"""
import numpy as np
import numexpr as ne
from scipy.special import roots_hermite
from ..base.constants import (
    InvbinSizePIhalf,
    binSize2,
)
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
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
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
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
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
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
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
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
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
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
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

def cross_section_PseudoVoigt(wn_grid, v, alpha, gamma, eta, coef, cutoff):
    """
    Calculate cross sections using Pseudo-Voigt line profile.

    Pseudo-Voigt is a weighted sum of Gaussian and Lorentzian profiles.

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
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            PseudoVoigt = PseudoVoigt_profile(dv, alpha, gamma, eta)
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
                _alpha = alpha[filter]
                _gamma = gamma[filter]
                _eta = eta[filter]
                _coef = coef[filter]
                PseudoVoigt = PseudoVoigt_profile(_dv, _alpha, _gamma, _eta)
                _xsec[idx] = _coef @ PseudoVoigt
    xsec[start:end] += _xsec   
    return xsec

def cross_section_BinnedGaussian(wn_grid, v, alpha, coef, cutoff):
    """
    Calculate cross sections using binned Gaussian line profile.

    Integrates Gaussian profile over bin width for each grid point.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array (bin centers), shape (n_points,)
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
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            BinnedGaussian = BinnedGaussian_profile(dv, alpha)
            _xsec[idx] = coef @ BinnedGaussian 
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
                BinnedGaussian = BinnedGaussian_profile(_dv, _alpha)
                _xsec[idx] = _coef @ BinnedGaussian 
    xsec[start:end] += _xsec / binSize2
    return xsec

def cross_section_BinnedLorentzian(wn_grid, v, gamma, coef, cutoff):
    """
    Calculate cross sections using binned Lorentzian line profile.

    Integrates Lorentzian profile over bin width for each grid point.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array (bin centers), shape (n_points,)
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
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
        wngrid_start = wn_grid[start]
        wngrid_end = wn_grid[end-1]
        bnormBinsize = ne.evaluate('1/(arctan((wngrid_end-v)/gamma)-arctan((wngrid_start-v)/gamma))/bin_size')
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            BinnedLorentzian = BinnedLorentzian_profile(dv, gamma, bnormBinsize)
            _xsec[idx] = coef @ BinnedLorentzian
    else:
        start = max(0,wn_grid.searchsorted(min(v)-cutoff)-1)
        end = min(wn_grid.searchsorted(max(v)+cutoff),len(wn_grid))
        wngrid_start = wn_grid[start]
        wngrid_end = wn_grid[end-1]
        bnormBinsize = ne.evaluate('1/(arctan((wngrid_end-v)/gamma)-arctan((wngrid_start-v)/gamma))/bin_size')
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            filter = np.abs(dv) <= cutoff
            if filter.sum() > 0:
                _dv = dv[filter]
                _gamma = gamma[filter]
                _bnormBinsize = bnormBinsize[filter]
                _coef = coef[filter]
                BinnedLorentzian = BinnedLorentzian_profile(_dv, _gamma, _bnormBinsize)
                _xsec[idx] = _coef @ BinnedLorentzian
    xsec[start:end] += _xsec  
    return xsec

def cross_section_BinnedVoigt(wn_grid, v, sigma, gamma, coef, cutoff):
    """
    Calculate cross sections using binned Voigt line profile.

    Uses Hermite quadrature to integrate Voigt profile over bin width.

    Parameters
    ----------
    wn_grid : np.ndarray
        Wavenumber grid array (bin centers), shape (n_points,)
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
    nquad = 20
    roots, weights = roots_hermite(nquad, mu=False)
    xsec = np.zeros_like(wn_grid)
    if cutoff == 'None':
        start = max(0,wn_grid.searchsorted(min(v))-1)
        end = min(wn_grid.searchsorted(max(v)),len(wn_grid))
        wngrid_start = wn_grid[start]
        wngrid_end = wn_grid[end-1]
        bnormqT = np.transpose(np.array([BinnedVoigt_bnormq(wngrid_start, wngrid_end, v, sigma, gamma, root) for root in roots]))
        _xsec = np.zeros(shape=(end-start))
        for i in range(start,end):
            idx = i-start
            dv = wn_grid[i] - v
            lorenzT = np.transpose(np.array([BinnedVoigt_lorenz(dv, sigma, gamma, root) for root in roots]))
            BinnedVoigtProfile = BinnedVoigt_profile(weights,bnormqT,lorenzT)                 
            _xsec[idx] = ne.evaluate('sum(coef * BinnedVoigtProfile)')
    else:
        start = max(0,wn_grid.searchsorted(min(v)-cutoff)-1)
        end = min(wn_grid.searchsorted(max(v)+cutoff),len(wn_grid))
        wngrid_start = wn_grid[start]
        wngrid_end = wn_grid[end-1]
        bnormq = [BinnedVoigt_bnormq(wngrid_start, wngrid_end, v, sigma, gamma, root) for root in roots]
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
                bnormqT = np.transpose(np.array([bnormq[i][filter] for i in range(nquad)]))
                lorenzT = np.transpose(np.array([BinnedVoigt_lorenz(_dv, _sigma, _gamma, root) for root in roots]))
                BinnedVoigtProfile = BinnedVoigt_profile(weights,bnormqT,lorenzT)    
                _xsec[idx] = ne.evaluate('sum(_coef * BinnedVoigtProfile)')
    xsec[start:end] += _xsec * InvbinSizePIhalf
    return xsec
