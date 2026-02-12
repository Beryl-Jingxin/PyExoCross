"""
Save modules for PyExoCross.

This package provides functions for saving calculation results to files.
Organized by database type (exomol, hitran).
"""
# ExoMol save functions
from .exomol.exomol_partition_func import save_exomol_partition_func
from .exomol.exomol_specific_heat import save_exomol_specific_heat
from .exomol.exomol_lifetime import save_exomol_lifetime
from .exomol.exomol_cooling_func import save_exomol_cooling_func
from .exomol.exomol_oscillator_strength import save_exomol_oscillator_strength
from .exomol.exomol_stick_spectra import save_exomol_stick_spectra
from .exomol.exomol_cross_section import save_exomol_cross_section
from .exomol.exomol_stick_spectra_cross_section import save_exomol_stick_spectra_cross_section

# HITRAN save functions
from .hitran.hitran_partition_func import save_hitran_partition_func
from .hitran.hitran_specific_heat import save_hitran_specific_heat
from .hitran.hitran_lifetime import save_hitran_lifetime
from .hitran.hitran_cooling_func import save_hitran_cooling_func
from .hitran.hitran_oscillator_strength import save_hitran_oscillator_strength
from .hitran.hitran_stick_spectra import save_hitran_stick_spectra
from .hitran.hitran_cross_section import save_hitran_cross_section
from .hitran.hitran_stick_spectra_cross_section import save_hitran_stick_spectra_cross_section

__all__ = [
    # ExoMol
    'save_exomol_partition_func',
    'save_exomol_specific_heat',
    'save_exomol_lifetime',
    'save_exomol_cooling_func',
    'save_exomol_oscillator_strength',
    'save_exomol_stick_spectra',
    'save_exomol_cross_section',
    'save_exomol_stick_spectra_cross_section',
    # HITRAN
    'save_hitran_partition_func',
    'save_hitran_specific_heat',
    'save_hitran_lifetime',
    'save_hitran_cooling_func',
    'save_hitran_oscillator_strength',
    'save_hitran_stick_spectra',
    'save_hitran_cross_section',
    'save_hitran_stick_spectra_cross_section',
]

