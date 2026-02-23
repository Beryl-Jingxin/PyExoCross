"""
Calculation modules for PyExoCross.

This package provides core calculation functions for:
- Absorption and emission coefficients
- Partition functions
- Specific heats
- Cooling functions
- Oscillator strengths
- Lifetimes
- Line profiles
- Parameter calculations
"""
from .calculate_intensity import (
    cal_abscoefs,
    cal_abscoefs_nlte_2T,
    cal_abscoefs_nlte_nvib,
    cal_abscoefs_nlte_pop,
)
from .calculate_emissivity import (
    cal_emicoefs,
    cal_emicoefs_nlte_2T,
    cal_emicoefs_nlte_nvib,
    cal_emicoefs_nlte_pop,
)
from .calculate_partition_func import (
    cal_partition_func,
    cal_Q_nlte_2T,
    cal_Q_nlte_nvib,
)
from .calculate_specific_heat import cal_specific_heat
from .calculate_cooling_func import cal_cooling_func
from .calculate_oscillator_strength import cal_oscillator_strength
from .calculate_lifetime import cal_lifetime
from .calculate_para import (
    cal_v,
    cal_Ep,
    cal_Ep_hitran,
    cal_Jp,
    cal_F,
    cal_uncertainty,
)
from .calcualte_line_profile import (
    Doppler_HWHM,
    Gaussian_standard_deviation,
    Lorentzian_HWHM,
    lifetime_broadening,
    DopplerHWHM_alpha,
    LorentzianHWHM_gamma,
    Doppler_profile,
    Lorentzian_profile,
    SciPyVoigt_profile,
    SciPyWofzVoigt_profile,
    Humlicek1,
    Humlicek2,
    Humlicek3,
    Humlicek4,
    HumlicekVoigt_profile,
    PseudoVoigt_profile,
    PseudoThompsonVoigt,
    PseudoKielkopfVoigt,
    PseudoOliveroVoigt,
    PseudoLiuLinVoigt,
    PseudoRoccoVoigt,
    BinnedGaussian_profile,
    BinnedLorentzian_profile,
    BinnedVoigt_bnormq,
    BinnedVoigt_lorenz,
    BinnedVoigt_profile,
)

__all__ = [
    # Absorption coefficients
    'cal_abscoefs',
    'cal_abscoefs_nlte_2T',
    'cal_abscoefs_nlte_nvib',
    'cal_abscoefs_nlte_pop',
    # Emission coefficients
    'cal_emicoefs',
    'cal_emicoefs_nlte_2T',
    'cal_emicoefs_nlte_nvib',
    'cal_emicoefs_nlte_pop',
    # Partition functions
    'cal_partition_func',
    'cal_Q_nlte_2T',
    'cal_Q_nlte_nvib',
    # Thermodynamic properties
    'cal_specific_heat',
    'cal_cooling_func',
    # Spectroscopic properties
    'cal_oscillator_strength',
    'cal_lifetime',
    # Parameter calculations
    'cal_v',
    'cal_Ep',
    'cal_Ep_hitran',
    'cal_Jp',
    'cal_F',
    'cal_uncertainty',
    # Line profiles
    'Doppler_HWHM',
    'Gaussian_standard_deviation',
    'Lorentzian_HWHM',
    'lifetime_broadening',
    'DopplerHWHM_alpha',
    'LorentzianHWHM_gamma',
    'Doppler_profile',
    'Lorentzian_profile',
    'SciPyVoigt_profile',
    'SciPyWofzVoigt_profile',
    'Humlicek1',
    'Humlicek2',
    'Humlicek3',
    'Humlicek4',
    'HumlicekVoigt_profile',
    'PseudoVoigt_profile',
    'PseudoThompsonVoigt',
    'PseudoKielkopfVoigt',
    'PseudoOliveroVoigt',
    'PseudoLiuLinVoigt',
    'PseudoRoccoVoigt',
    'BinnedGaussian_profile',
    'BinnedLorentzian_profile',
    'BinnedVoigt_bnormq',
    'BinnedVoigt_lorenz',
    'BinnedVoigt_profile',
]

