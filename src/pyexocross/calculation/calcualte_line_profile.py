"""
Calculate various line profiles for spectroscopic calculations.

This module provides functions for calculating Doppler, Lorentzian, Voigt,
and other line profiles used in cross-section calculations.
"""
import numpy as np
import numexpr as ne
from scipy.special import voigt_profile, wofz, erf, roots_hermite
from ..base.constants import (
    get_doppler_constants,
    Tref,
    Pref,
    c,
    PI,
    PI4c,
    Negln2,
    SqrtPI,
    Sqrtln2,
    InvSqrt2,
    InvSqrtPi,
    InvSqrt2Pi,
    InvSprtln2,
    InvSqrt2ln2,
    TwoSqrt2ln2,
    Sqrtln2InvPi,
    Sqrt2NAkBln2mInvc,
    binSizeHalf
)
from pyexocross.base.config_manager import get_config

# Lazy config: only resolved when line-profile functions run (avoids requiring
# config in worker processes that only import this module for other symbols).
_line_profile_config = None


def _get_line_profile_config():
    global _line_profile_config
    if _line_profile_config is None:
        _line_profile_config = get_config()
    return _line_profile_config

def Doppler_HWHM(v, T, mass=None):
    """
    Calculate the Doppler half-width at half-maximum (HWHM).

    alpha = sqrt(2 * N_A * kB * T * ln(2) / mass) * v / c

    Parameters
    ----------
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    T : float or np.ndarray
        Temperature in Kelvin
    mass : float
        Molecular or atomic mass in unit (Da)

    Returns
    -------
    np.ndarray
        Doppler HWHM (alpha) array, shape (n_levels,)
    """
    # If no explicit mass is provided, use the precomputed module constant
    # from pyexocross.base.constants (legacy behavior). Otherwise compute from mass.
    if mass is None:
        doppler_const = Sqrt2NAkBln2mInvc
    else:
        doppler_const = get_doppler_constants(mass)
    alpha = ne.evaluate('doppler_const * sqrt(T) * v')
    return alpha

def Gaussian_standard_deviation(alpha):
    """
    Calculate the Gaussian standard deviation from HWHM.

    sigma = alpha / sqrt(2 * ln(2))

    Parameters
    ----------
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Gaussian standard deviation (sigma) array, shape (n_levels,)
    """
    sigma = ne.evaluate('alpha * InvSqrt2ln2')
    return sigma

def Lorentzian_HWHM(gamma_L, n_air,T,P):
    """
    Calculate the Lorentzian half-width at half-maximum (HWHM).

    gamma = gamma_L * (Tref / T)^n_air * (P / Pref)

    Parameters
    ----------
    gamma_L : np.ndarray
        Reference Lorentzian HWHM array, shape (n_levels,)
    n_air : np.ndarray
        Temperature exponent array, shape (n_levels,)
    T : float or np.ndarray
        Temperature in Kelvin
    P : float or np.ndarray
        Pressure in bar

    Returns
    -------
    np.ndarray
        Lorentzian HWHM (gamma) array, shape (n_levels,)
    """
    # Ensure T and P are numpy arrays/scalars that can broadcast correctly
    T_arr = np.asarray(T)
    P_arr = np.asarray(P)
    gamma = ne.evaluate('gamma_L * (Tref / T_arr)**n_air * (P_arr / Pref)')
    return gamma

def lifetime_broadening(tau):
    """
    Calculate the lifetime broadening contribution to linewidth.

    gamma_tau = 1 / (4 * pi * c) * tau

    Parameters
    ----------
    tau : np.ndarray
        Lifetime array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Lifetime broadening (gamma_tau) array, shape (n_levels,)
    """
    gamma_tau = ne.evaluate('1 / (PI4c) * tau')
    return gamma_tau

def DopplerHWHM_alpha(num_v, alpha_HWHM, v, T):
    """
    Calculate Doppler HWHM based on user settings.

    Parameters
    ----------
    num_v : int
        Number of transitions
    alpha_HWHM : float or np.ndarray
        User-provided Doppler HWHM value or array
    v : np.ndarray
        Wavenumber array, shape (n_levels,)
    T : float or np.ndarray
        Temperature in Kelvin

    Returns
    -------
    np.ndarray
        Doppler HWHM array, shape (n_levels,)
    """
    (DopplerHWHMYN,) = _get_line_profile_config()[:1]
    if DopplerHWHMYN == 'Y' and num_v > 0:
        alpha = np.full(num_v, alpha_HWHM)
    elif DopplerHWHMYN == 'U':
        alpha = alpha_HWHM
    else:
        alpha = Doppler_HWHM(v,T)
    return alpha

def LorentzianHWHM_gamma(num, gamma_HWHM, nbroad, gamma_L, n_air, gamma_air, gamma_self, tau, T, P):
    """
    Calculate Lorentzian HWHM based on database and user settings.

    Parameters
    ----------
    num : int
        Number of transitions
    gamma_HWHM : float or np.ndarray
        User-provided Lorentzian HWHM value or array
    nbroad : int
        Number of broadening species
    gamma_L : dict or np.ndarray
        Reference Lorentzian HWHM (dict for ExoMol, array for HITRAN)
    n_air : dict or np.ndarray
        Temperature exponent (dict for ExoMol, array for HITRAN)
    gamma_air : np.ndarray
        Air-broadening coefficient (HITRAN only)
    gamma_self : np.ndarray
        Self-broadening coefficient (HITRAN only)
    tau : np.ndarray
        Lifetime array for predissociation broadening
    T : float or np.ndarray
        Temperature in Kelvin
    P : float or np.ndarray
        Pressure in bar

    Returns
    -------
    np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Raises
    ------
    ValueError
        If gamma_L or n_air array lengths don't match num
    """
    (
        _,
        LorentzianHWHMYN,
        database,
        predissocYN,
        check_predissoc,
        profile,
    ) = _get_line_profile_config()
    if LorentzianHWHMYN == 'Y' and num > 0:
        gamma = np.full(num, gamma_HWHM)
    elif LorentzianHWHMYN == 'U':
        gamma = gamma_HWHM
    elif database == 'ExoMol' or database == 'ExoAtom' and num > 0:
        # gamma_L and n_air are now dictionaries with numpy arrays
        gamma_contributions = []
        for i in range(nbroad):
            # Get arrays directly from dictionary
            gamma_L_i = np.asarray(gamma_L[i])
            n_air_i = np.asarray(n_air[i])
            
            # Ensure they are 1D arrays of correct length
            gamma_L_i = gamma_L_i.flatten() if gamma_L_i.ndim > 1 else gamma_L_i
            n_air_i = n_air_i.flatten() if n_air_i.ndim > 1 else n_air_i
            
            # Verify shapes match
            if len(gamma_L_i) != num:
                raise ValueError(f"gamma_L[{i}] has length {len(gamma_L_i)}, expected {num}")
            if len(n_air_i) != num:
                raise ValueError(f"n_air[{i}] has length {len(n_air_i)}, expected {num}")
            
            gamma_contributions.append(Lorentzian_HWHM(gamma_L_i, n_air_i, T, P))
        gamma = sum(gamma_contributions)
        if predissocYN == 'Y' and check_predissoc == 0 and 'VOI' in profile:
            gamma += lifetime_broadening(tau)
        else:
            pass
    elif database == 'HITRAN' or database == 'HITEMP' and num > 0:  
        gamma_L = gamma_air*0.7 + gamma_self*0.3
        gamma = Lorentzian_HWHM(gamma_L, n_air,T,P) 
    else:
        gamma = np.zeros(num)
    return gamma

# def FWHM(alpha, gamma):
#     '''Return the Gaussian full-width at half-maximum (FWHM)   -- fG.
#        Return the Lorentzian full-width at half-maximum (FWHM) -- fL.
#     '''
#     # fG = 2 * sigma * np.sqrt(2 * np.log(2)) = 2 * alpha
#     # fG = ne.evaluate('sigma * TwoSqrt2ln2')
#     fG = ne.evaluate('2 * alpha')
#     fL = ne.evaluate('2 * gamma')
#     return (fG, fL)

def Doppler_profile(dv, alpha):
    """
    Calculate Doppler (Gaussian) line profile.

    Doppler = sqrt(ln(2) / pi) / alpha * exp(-ln(2) * (dv / alpha)^2) 

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center, shape (n_points,)
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Doppler line profile array, shape (n_levels, n_points)
    """
    DopplerProfile = ne.evaluate('Sqrtln2InvPi / alpha * exp(Negln2 * (dv / alpha)**2)')
    return DopplerProfile

def Lorentzian_profile(dv, gamma):
    """
    Calculate Lorentzian line profile.

    Lorentzian = gamma / pi / (dv^2 + gamma^2)

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center, shape (n_points,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Lorentzian line profile array, shape (n_levels, n_points)
    """
    LorentzianProfile = ne.evaluate('gamma / PI / (dv**2 + gamma**2)')
    return LorentzianProfile

def SciPyVoigt_profile(dv, sigma, gamma):
    """
    Calculate Voigt line profile using scipy.special.voigt_profile.

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center, shape (n_points,)
    sigma : np.ndarray
        Gaussian standard deviation array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Voigt line profile array, shape (n_levels, n_points)
    """
    SciPyVoigtProfile = voigt_profile(dv, sigma, gamma)
    return SciPyVoigtProfile

def SciPyWofzVoigt_profile(dv, sigma, gamma):
    """
    Calculate Voigt line profile using scipy.special.wofz (Faddeeva function).

    SciPyWofzVoigt = real(wofz((dv + 1j*gamma)/sigma/sqrt(2))) / sigma / sqrt(2*pi)

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center, shape (n_points,)
    sigma : np.ndarray
        Gaussian standard deviation array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Voigt line profile array, shape (n_levels, n_points)
    """
    z = ne.evaluate('(dv + 1j*gamma)/sigma*InvSqrt2')
    wz = wofz(z)
    SciPyWofzVoigtProfile = ne.evaluate('real(wz) / sigma * InvSqrt2Pi')
    return SciPyWofzVoigtProfile

def Humlicek1(t):
    """
    Humlicek rational approximation region 1 for Voigt profile.

    w = t / (0.5 + t^2) / sqrt(pi) 

    Parameters
    ----------
    t : np.ndarray
        Complex argument t = y - i*x

    Returns
    -------
    np.ndarray
        Complex Faddeeva function approximation
    """
    w = ne.evaluate('t * InvSqrtPi / (0.5 + t**2)')
    return(w)

def Humlicek2(t, u):
    """
    Humlicek rational approximation region 2 for Voigt profile.

    Parameters
    ----------
    t : np.ndarray
        Complex argument t = y - i*x
    u : np.ndarray
        Complex argument u = t^2

    Returns
    -------
    np.ndarray
        Complex Faddeeva function approximation
    """
    w = ne.evaluate('(t*(1.4104739589+u*InvSqrtPi))/(0.75+u*(3+u))')
    return(w)

def Humlicek3(t):
    """
    Humlicek rational approximation region 3 for Voigt profile.

    Parameters
    ----------
    t : np.ndarray
        Complex argument t = y - i*x

    Returns
    -------
    np.ndarray
        Complex Faddeeva function approximation
    """
    w = ne.evaluate('(16.4955+t*(20.20933+t*(11.96482+t*(3.778987+0.5642236*t))))/(16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))))')
    return(w)

def Humlicek4(t, u):
    """
    Humlicek rational approximation region 4 for Voigt profile.

    Parameters
    ----------
    t : np.ndarray
        Complex argument t = y - i*x
    u : np.ndarray
        Complex argument u = t^2

    Returns
    -------
    np.ndarray
        Complex Faddeeva function approximation
    """
    nom = ne.evaluate('t*(36183.31-u*(3321.99-u*(1540.787-u*(219.031-u*(35.7668-u*(1.320522-u*0.56419))))))')
    den = ne.evaluate('32066.6-u*(24322.8-u*(9022.23-u*(2186.18-u*(364.219-u*(61.5704-u*(1.84144-u))))))')
    w = ne.evaluate('exp(u)-nom/den')   
    return(w)

def HumlicekVoigt_profile(dv, alpha, gamma):
    """
    Calculate Voigt line profile using Humlicek's rational approximation.

    Uses four-region approximation for efficient computation of Faddeeva function.

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center, shape (n_points,)
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Voigt line profile array, shape (n_levels, n_points)
    """
    x = ne.evaluate('dv * Sqrtln2 / alpha')
    y = ne.evaluate('gamma * Sqrtln2 / alpha')
    t = ne.evaluate('y-1j*x')
    s = ne.evaluate('abs(x)+y')
    u = ne.evaluate('t**2')
    w = np.zeros_like(s)
    ybound = ne.evaluate('0.195*abs(x)-0.176')
    # Region 1
    humfilter1 = s >= 15
    t1 = t[humfilter1]
    w[humfilter1] = Humlicek1(t1)
    # Region 2
    humfilter2 = (5.5 <= s) & (s < 15)
    t2 = t[humfilter2]
    u2 = u[humfilter2]
    w[humfilter2] = Humlicek2(t2, u2)
    # Region 3
    humfilter3 = (s < 5.5) & (y >= ybound)
    t3 = t[humfilter3]
    w[humfilter3] = Humlicek3(t3)
    # Region 4
    humfilter4 = (s < 5.5) & (y < ybound)
    t4 = t[humfilter4]
    u4 = u[humfilter4]
    w[humfilter4] = Humlicek4(t4, u4)  
    HumlicekVoigtProfile = ne.evaluate('real(w) / alpha * Sqrtln2InvPi')
    return HumlicekVoigtProfile

def PseudoVoigt_profile(dv, alpha, gamma, eta):
    """
    Calculate Pseudo-Voigt line profile as weighted sum of Gaussian and Lorentzian.

    PseudoVoigt = eta * Lorentzian + (1 - eta) * Gaussian

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center, shape (n_points,)
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)
    eta : np.ndarray
        Mixing parameter array (0 = pure Gaussian, 1 = pure Lorentzian), shape (n_levels,)

    Returns
    -------
    np.ndarray
        Pseudo-Voigt line profile array, shape (n_levels, n_points)
    """
    GaussianProfile = Doppler_profile(dv, alpha)
    LorentzianProfile = Lorentzian_profile(dv, gamma)
    PseudoVoigtProfile = ne.evaluate('eta * LorentzianProfile + (1 - eta) * GaussianProfile')
    return PseudoVoigtProfile

def PseudoThompsonVoigt(alpha, gamma):
    """
    Calculate mixing parameter eta for Pseudo-Voigt using Thompson's approximation.

    hv is the Voigt half-width at half-maximum (HWHM), which can be found from the HWHM of the associated Doppler and Lorentzian profile.

    Parameters
    ----------
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        eta array is a function of Voigt profile half width at half maximum (HWHM) parameters, shape (n_levels,)
    """ 
    hV = ne.evaluate('(alpha**5+2.69269*alpha**4*gamma+2.42843*alpha**3*gamma**2+4.47163*alpha**2*gamma**3+0.07842*alpha*gamma**4+gamma**5)**0.2')
    eta = ne.evaluate('1.36603*(gamma/hV) - 0.47719*(gamma/hV)**2 + 0.11116*(gamma/hV)**3')
    return (eta)

def PseudoKielkopfVoigt(alpha, gamma):
    """
    Calculate mixing parameter eta for Pseudo-Voigt using Kielkopf's approximation.

    hv is the Voigt half-width at half-maximum (HWHM), which can be found from the HWHM of the associated Doppler and Lorentzian profile.

    Parameters
    ----------
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        eta array is a function of Voigt profile half width at half maximum (HWHM) parameters, shape (n_levels,)
    """
    hV = ne.evaluate('0.5346 * gamma + sqrt(0.2166 * gamma**2 + alpha**2)')
    eta = ne.evaluate('1.36603*(gamma/hV) - 0.47719*(gamma/hV)**2 + 0.11116*(gamma/hV)**3')
    return (eta)

def PseudoOliveroVoigt(alpha, gamma):
    """
    Calculate mixing parameter eta for Pseudo-Voigt using Olivero's approximation.

    hv is the Voigt half-width at half-maximum (HWHM), which can be found from the HWHM of the associated Doppler and Lorentzian profile.

    Parameters
    ----------
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        eta array is a function of Voigt profile half width at half maximum (HWHM) parameters, shape (n_levels,)
    """
    d = ne.evaluate('(gamma-alpha)/(gamma+alpha)')
    hV = ne.evaluate('(1-0.18121*(1-d**2)-(0.023665*exp(0.6*d)+0.00418*exp(-1.9*d))*sin(PI*d))*(alpha+gamma)')
    eta = ne.evaluate('1.36603*(gamma/hV) - 0.47719*(gamma/hV)**2 + 0.11116*(gamma/hV)**3')
    return (eta)

def PseudoLiuLinVoigt(alpha, gamma):
    """
    Calculate mixing parameter eta for Pseudo-Voigt using Liu-Lin's approximation.

    hv is the Voigt half-width at half-maximum (HWHM), which can be found from the HWHM of the associated Doppler and Lorentzian profile.

    Parameters
    ----------
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        eta array is a function of Voigt profile half width at half maximum (HWHM) parameters, shape (n_levels,)
    """
    d = ne.evaluate('(gamma-alpha)/(gamma+alpha)')
    hV = ne.evaluate('(1-0.18121*(1-d**2)-(0.023665*exp(0.6*d)+0.00418*exp(-1.9*d))*sinPI*d)*(alpha+gamma)')
    eta = ne.evaluate('0.68188+0.61293*d-0.18384*d**2-0.11568*d**3')
    return (eta)

def PseudoRoccoVoigt(alpha, gamma):
    """
    Calculate mixing parameter eta for Pseudo-Voigt using Rocco's approximation.

    hv is the Voigt half-width at half-maximum (HWHM), which can be found from the HWHM of the associated Doppler and Lorentzian profile.

    Parameters
    ----------
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        eta array is a function of Voigt profile half width at half maximum (HWHM) parameters, shape (n_levels,)
    """
    y = gamma*Sqrtln2/alpha
    erfy = erf(y)
    bhalfy = ne.evaluate('y+Sqrtln2*exp(-0.6055*y+0.0718*y**2-0.0049*y**3+0.000136*y**4)')
    Vy = ne.evaluate('bhalfy*exp(y**2)*(1-erfy)')
    hV = ne.evaluate('alpha / Sqrtln2 * bhalfy')
    eta = ne.evaluate('(Vy-Sqrtln2)/(Vy*OneminSqrtPIln2)')
    return (eta)

def BinnedGaussian_profile(dv, alpha):
    """
    Calculate binned Gaussian line profile integrated over bin width.

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center (bin center), shape (n_points,)
    alpha : np.ndarray
        Doppler HWHM array, shape (n_levels,)

    Returns
    -------
    np.ndarray
        Binned Gaussian line profile array, shape (n_levels, n_points)
    """
    erfxpos = erf(ne.evaluate('Sqrtln2*(dv+binSizeHalf)/alpha'))
    erfxneg = erf(ne.evaluate('Sqrtln2*(dv-binSizeHalf)/alpha'))
    BinnedGaussianProfile = ne.evaluate('erfxpos-erfxneg')
    return BinnedGaussianProfile

def BinnedLorentzian_profile(dv, gamma, bnormBinsize):
    """
    Calculate binned Lorentzian line profile integrated over bin width.

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center (bin center), shape (n_points,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)
    bnormBinsize : float
        Normalization factor for bin size

    Returns
    -------
    np.ndarray
        Binned Lorentzian line profile array, shape (n_levels, n_points)
    """
    BinnedLorentzianProfile = ne.evaluate('(arctan((dv+binSizeHalf)/gamma)-arctan((dv-binSizeHalf)/gamma))*bnormBinsize')
    return BinnedLorentzianProfile

def BinnedVoigt_bnormq(wngrid_start, wngrid_end, v, sigma, gamma, x):
    """
    Calculate normalization factor for binned Voigt profile.

    Parameters
    ----------
    wngrid_start : float
        Start of wavenumber grid bin
    wngrid_end : float
        End of wavenumber grid bin
    v : np.ndarray
        Line center wavenumber array, shape (n_levels,)
    sigma : np.ndarray
        Gaussian standard deviation array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)
    x : np.ndarray
        Hermite quadrature points, shape (n_points,)

    Returns
    -------
    np.ndarray
        Normalization factor array, shape (n_levels, n_points)
    """
    vxsigma = ne.evaluate('v+x*sigma')
    bnormq = ne.evaluate('1/(arctan((wngrid_end-vxsigma)/gamma)-arctan((wngrid_start-vxsigma)/gamma))')
    return bnormq

def BinnedVoigt_lorenz(dv, sigma, gamma, x):
    """
    Calculate Lorentzian component for binned Voigt profile.

    Parameters
    ----------
    dv : np.ndarray
        Frequency offset from line center, shape (n_points,)
    sigma : np.ndarray
        Gaussian standard deviation array, shape (n_levels,)
    gamma : np.ndarray
        Lorentzian HWHM array, shape (n_levels,)
    x : np.ndarray
        Hermite quadrature points, shape (n_quad_points,)

    Returns
    -------
    np.ndarray
        Lorentzian component array, shape (n_levels, n_points, n_quad_points)
    """
    dvxsigma = ne.evaluate('dv-x*sigma')
    lorenz = ne.evaluate('arctan((dvxsigma+binSizeHalf)/gamma)-arctan((dvxsigma-binSizeHalf)/gamma)')
    return lorenz

def BinnedVoigt_profile(w, bnormq, lorenz):
    """
    Calculate binned Voigt line profile using Hermite quadrature.

    Parameters
    ----------
    w : np.ndarray
        Hermite quadrature weights, shape (n_quad_points,)
    bnormq : np.ndarray
        Normalization factor array, shape (n_levels, n_quad_points)
    lorenz : np.ndarray
        Lorentzian component array, shape (n_levels, n_points, n_quad_points)

    Returns
    -------
    np.ndarray
        Binned Voigt line profile array, shape (n_levels, n_points)
    """
    BinnedVoigtProfile = ne.evaluate('sum(w*bnormq*lorenz,axis=1)')
    return BinnedVoigtProfile

# Line Profile
def line_profile(profile):
    """
    Convert profile code string to human-readable profile label.

    Maps profile codes (e.g., 'DOP', 'GAU', 'LOR', 'SCI', etc.) to
    descriptive profile names.

    Parameters
    ----------
    profile : str
        Profile code string

    Returns
    -------
    str
        Human-readable profile label (e.g., 'Doppler', 'Gaussian', 'Lorentzian', 'SciPy Voigt')

    Raises
    ------
    ValueError
        If profile code is not recognized
    """
    if profile[0:3] == 'DOP':
        profile_label = 'Doppler'
    elif profile[0:3] == 'GAU':
        profile_label = 'Gaussian'
    elif profile[0:3] == 'LOR':
        profile_label = 'Lorentzian'
    elif 'SCI' in profile and 'W' not in profile:
        profile_label = 'SciPy Voigt'
    elif 'W' in profile:
        profile_label = 'SciPy wofz Voigt'
    elif 'H' in profile:
        profile_label = 'Humlicek Voigt'  
    elif 'TH' in profile:
        profile_label = 'Thompson pseudo-Voigt'
    elif 'K' in profile and 'H' not in profile:
        profile_label = 'Kielkopf pseudo-Voigt'
    elif 'OL' in profile:
        profile_label = 'Olivero pseudo-Voigt'
    elif 'LI' in profile or 'LL' in profile:
        profile_label = 'Liu-Lin pseudo-Voigt'
    elif 'RO' in profile:
        profile_label = 'Rocco pseudo-Voigt'
    elif 'BIN' in profile and 'DOP' in profile:
        profile_label = 'Binned Doppler'     
    elif 'BIN' in profile and 'GAU' in profile:
        profile_label = 'Binned Gaussion'
    elif 'BIN' in profile and 'LOR' in profile:
        profile_label = 'Binned Lorentzian'
    elif 'BIN' in profile and 'VOI' in profile:
        profile_label = 'Binned Voigt'
    else:
        raise ValueError('Please choose line profile from the list.')
    return profile_label
