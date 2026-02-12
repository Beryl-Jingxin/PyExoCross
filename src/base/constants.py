"""
Physical constants and calculation parameters for PyExoCross.

This module provides all physical constants and derived constants used throughout
the codebase. Some constants are computed dynamically based on molecular mass
and bin size parameters.
"""
import numpy as np
import multiprocessing as mp
import astropy.constants as ac

# File processing constants
DEFAULT_CHUNK_SIZE = 100_000  # fallback before user input is parsed
LARGE_TRANS_FILE_BYTES = 1024 ** 3  # 1 GB threshold for large .trans files
MAX_LARGE_FILE_WORKERS = 2
MAX_INFLIGHT_MULTIPLIER = 2
LARGE_WRITE_CHUNK_ROWS = 200_000
num_cpus = mp.cpu_count()

# Physical constants
Tref = 296.0                        # Reference temperature is 296 K
Pref = 1.0                          # Reference pressure is 1 bar
N_A = ac.N_A.value                  # Avogadro number (1/mol)
h = ac.h.to('erg s').value          # Planck's const (erg s)
c = ac.c.to('cm/s').value           # Velocity of light (cm/s)
kB = ac.k_B.to('erg/K').value       # Boltzmann's const (erg/K)
R = ac.R.to('J / (K mol)').value    # Molar gas constant (J/(K mol))
c2 = h * c / kB                     # Second radiation constant (cm K)

c2InvTref = c2 / Tref               # c2 / T_ref (cm)    # erg cm
PI = np.pi
hc = h * c
ln22 = np.log(2) * 2
sinPI = np.sin(np.pi)
SqrtPI = np.sqrt(np.pi)
Sqrtln2 = np.sqrt(np.log(2))
OneminSqrtPIln2 = 1 - np.sqrt(np.pi * np.log(2))
Negln2 = -np.log(2)
PI4c = np.pi * 4 * c
Inv8Pic = 1 / (8 * np.pi * c)       # 1 / (8 * pi * c) (s/cm)
hcInv4Pi = hc / (4 * np.pi)
Inv2ln2 = 1 / (2 * np.log(2))
InvSqrt2 = 1 / np.sqrt(2)
InvSqrtPi = 1 / np.sqrt(np.pi)
InvSprtln2 = 1 / np.sqrt(np.log(2))
InvSqrt2Pi = 1 / np.sqrt(2 * np.pi)
InvSqrt2ln2 = 1 / np.sqrt(2 * np.log(2))
TwoSqrt2ln2 = 2 * np.sqrt(2 * np.log(2))
Sqrtln2InvPi = np.sqrt(np.log(2) / np.pi)

# Legacy/derived constants expected by older modules.
# They are initialized with safe defaults and updated later once bin_size
# and mass information are available.
#
# Doppler-related constant (updated per-mass inside Doppler_HWHM in
# calcualte_line_profile.py; this default is only to satisfy imports).
Sqrt2NAkBln2mInvc = 0.0

# Bin-size–related constants. These are updated from bin_size via
# get_bin_size_constants(...) inside the core execution flow.
binSize2 = 1.0
binSizeHalf = 0.5
InvbinSizePIhalf = 1.0 / (1.0 * np.sqrt(np.pi))


def get_doppler_constants(mass):
    """
    Calculate Doppler broadening constants for a given molecular mass.

    Parameters
    ----------
    mass : float
        Molecular mass in atomic mass units (Da)

    Returns
    -------
    float
        Doppler constant: sqrt(2 * N_A * kB * ln(2) / mass) / c
    """
    return np.sqrt(2 * N_A * kB * np.log(2) / mass) / c


def get_bin_size_constants(bin_size):
    """
    Calculate bin size related constants.

    Parameters
    ----------
    bin_size : float
        Bin size value

    Returns
    -------
    dict
        Dictionary containing bin size constants:
        - binSize2: bin_size * 2
        - binSizePI: bin_size * pi
        - binSizeHalf: bin_size / 2
        - InvbinSizePIhalf: 1 / (bin_size * sqrt(pi))
    """
    if bin_size == 'None' or bin_size is None:
        return {}
    return {
        'binSize2': bin_size * 2,
        'binSizePI': bin_size * np.pi,
        'binSizeHalf': bin_size / 2,
        'InvbinSizePIhalf': 1 / (bin_size * np.sqrt(np.pi))
    }