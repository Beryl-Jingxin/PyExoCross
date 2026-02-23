"""
Test all PyExoCross API functions using ExoMol MgH parameters.

Parameters are derived from: .input/MgH_ExoMol.inp
"""
import sys
import os

# Ensure src layout is importable in local runs
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src_root = os.path.join(project_root, 'src')
if src_root not in sys.path:
    sys.path.insert(0, src_root)

import pyexocross as px


# ---------------------------------------------------------------------------
# Common parameters (shared across all functions)
# ---------------------------------------------------------------------------
COMMON = dict(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/Users/beryl/Academic/UCL/PhD/Data/database/ExoAtom/', #'/home/jingxin/LHD/Program/Databases/ExoAtom/',
    save_path='/Users/beryl/Academic/UCL/PhD/Data/pyexocross/', #'/home/jingxin/LHD/Program/Data/pyexocross/',
    logs_path='/Users/beryl/Academic/UCL/PhD/Data/pyexocross/log/test_api_exoatom.log', #'/home/jingxin/LHD/Program/Data/pyexocross/log/test_api_exoatom.log', 
)

# Quantum number labels/formats (needed by conversion, stick_spectra, cross_sections)
QN_PARAMS = dict(
    qnslabel_list=['configuration', 'Multiple', 'parity'],
    qnsformat_list=['%20s', '%10s', '%2s'],
)

# NLTE parameters (needed by stick_spectra and cross_sections)
NLTE_PARAMS = dict(
    # Non-LTE: Population
    nlte_method='P',
    nlte_path='/home/jingxin/LHD/Program/Databases/ExoAtom/Ar/NIST/Ar_Ids.csv',
)

# Spectral range parameters
RANGE_PARAMS = dict(
    temperatures=[1000, 2000],      # Temperature in unit of K
    wn_wl='WN',                     # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
    wn_wl_unit='cm-1',              # Unit for wavenumber (default: cm⁻¹)
    min_range=0,                    # Minimum wavenumber in unit of cm⁻¹
    max_range=115400,               # Maximum wavenumber in unit of cm⁻¹
    abs_emi='Ab',                   # Absorption or emission (default: 'Absorption')
    unc_filter=None,                # Uncertainty filter (default: None)
    threshold=None,                 # Threshold filter (default: None)
)

# Cores and chunks
COMPUTE_PARAMS = dict(
    ncputrans=4,                    # Number of CPU threads for each transition file (default: 4)
    ncpufiles=1,                    # Number of CPU files for transition calculation (default: 1)
    chunk_size=10000,               # Chunk size for transition calculation (default: 100000)
)


# ---------------------------------------------------------------------------
# Test functions
# ---------------------------------------------------------------------------
def test_conversion():
    """Test ExoMol -> HITRAN conversion."""
    print('\n' + '='*70)
    print('TEST: px.conversion()  [ExoMol -> HITRAN]')
    print('='*70)
    px.conversion(
        **COMMON,
        **QN_PARAMS,
        **COMPUTE_PARAMS,
        conversion_format=1,
        conversion_min_freq=0,          # Minimum wavenumber in unit of cm⁻¹
        conversion_max_freq=115400,     # Maximum wavenumber in unit of cm⁻¹
        conversion_unc=None,            # Uncertainty filter (default: None)
        conversion_threshold=None,      # Threshold filter (default: None)
        global_qn_label_list=['configuration'],             # Quantum number label for global quantum numbers
        global_qn_format_list=['%20s'],                     # Quantum number format for global quantum numbers
        local_qn_label_list=['J', 'Multiple', 'parity'],    # Quantum number label for local quantum numbers
        local_qn_format_list=['%5.1f', '%10s', '%1s'],      # Quantum number format for local quantum numbers
    )
    print('PASSED: conversion()')


def test_partition_functions():
    """Test partition function calculation."""
    print('\n' + '='*70)
    print('TEST: px.partition_functions()')
    print('='*70)
    px.partition_functions(
        **COMMON,
        **COMPUTE_PARAMS,
        ntemp=1,                     # Number of temperature steps in unit of K (default: 1)
        tmax=6000,                   # Maximum temperature in unit of K (default: 5000)
    )
    print('PASSED: partition_functions()')


def test_specific_heats():
    """Test specific heat calculation."""
    print('\n' + '='*70)
    print('TEST: px.specific_heats()')
    print('='*70)
    px.specific_heats(
        **COMMON,
        **COMPUTE_PARAMS,
        ntemp=1,                     # Number of temperature steps in unit of K (default: 1)
        tmax=6000,                   # Maximum temperature in unit of K (default: 5000)
    )
    print('PASSED: specific_heats()')


def test_cooling_functions():
    """Test cooling function calculation."""
    print('\n' + '='*70)
    print('TEST: px.cooling_functions()')
    print('='*70)
    px.cooling_functions(
        **COMMON,
        **COMPUTE_PARAMS,
        ntemp=1,                     # Number of temperature steps in unit of K (default: 1)
        tmax=6000,                   # Maximum temperature in unit of K (default: 5000)
    )
    print('PASSED: cooling_functions()')


def test_lifetimes():
    """Test radiative lifetime calculation."""
    print('\n' + '='*70)
    print('TEST: px.lifetimes()')
    print('='*70)
    px.lifetimes(
        **COMMON,
        **COMPUTE_PARAMS,
        compress=False,               # Whether to compress the states file (default: False)
    )
    print('PASSED: lifetimes()')


def test_oscillator_strengths():
    """Test oscillator strength calculation."""
    print('\n' + '='*70)
    print('TEST: px.oscillator_strengths()')
    print('='*70)
    px.oscillator_strengths(
        **COMMON,
        **COMPUTE_PARAMS,
        gf_or_f='f',                  # 'gf' for weighted oscillator strength, 'f' for f-value (default: 'f')
        plot=True,                    # Whether to plot results (default: False)
        plot_method='log',            # Plot in linear (lin) or logarithm (log) (default: 'log')
        plot_wn_wl='WN',              # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='cm-1',             # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,            # Lower limit for y-axis (default: 1e-30)
    )
    print('PASSED: oscillator_strengths()')


def test_stick_spectra():
    """Test stick spectra calculation."""
    print('\n' + '='*70)
    print('TEST: px.stick_spectra()')
    print('='*70)
    px.stick_spectra(
        **COMMON,
        **QN_PARAMS,
        **NLTE_PARAMS,          # If Non-LTE is enabled, this parameter is required.
        **RANGE_PARAMS,
        **COMPUTE_PARAMS,
        qns_filter={
            'configuration': [],
            'Multiple': [],
            'parity': [],
        },
        plot=True,              # Whether to plot results (default: False)
        plot_method='log',      # Plot in linear (lin) or logarithm (log) (default: 'log')
        plot_wn_wl='Wn',        # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='cm-1',       # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,      # Lower limit for y-axis (default: 1e-30 cm/molecule)
    )
    print('PASSED: stick_spectra()')


def test_cross_sections():
    """Test cross section calculation."""
    print('\n' + '='*70)
    print('TEST: px.cross_sections()')
    print('='*70)
    px.cross_sections(
        **COMMON,
        **QN_PARAMS,
        **NLTE_PARAMS,          # If Non-LTE is enabled, this parameter is required.
        **RANGE_PARAMS,
        **COMPUTE_PARAMS,
        pressures=[1.0],        # Pressure in bar (default: [1.0])
        bin_size=0.1,           # Bin size for wavenumber grid
        profile='SciPyVoigt',   # Line profile name (default: 'Gaussian')
        predissociation=False,  # Predissociation (default: False)
        cutoff=None,            # Cutoff distance in cm⁻¹ (default: None)
        broadeners=['Default'], # Broadening species (default: ['Default'])
        ratios=[1.0],           # Broadening ratios (default: [1.0])
        alpha_hwhm=None,        # Constant Doppler HWHM (None, will calculate from broadening) or custom value (default: 3.0)
        gamma_hwhm=None,        # Constant Lorentzian HWHM (None, will calculate from broadening) or custom value (default: 0.5)
        plot=True,              # Whether to plot results (default: False)
        plot_method='log',      # Plot in linear (lin) or logarithm (log)
        plot_wn_wl='WN',        # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='cm-1',       # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,      # Lower limit for y-axis (default: 1e-30 cm²/molecule)
    )
    print('PASSED: cross_sections()')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print('PyExoCross API Test — ExoMol MgH')
    print(f'pyexocross version: {px.__version__}')

    tests = [
        test_conversion,
        test_partition_functions,
        test_specific_heats,
        test_cooling_functions,
        test_lifetimes,
        test_oscillator_strengths,
        test_stick_spectra,
        test_cross_sections,
    ]

    passed = 0
    failed = 0
    for test_fn in tests:
        try:
            test_fn()
            passed += 1
        except Exception as exc:
            failed += 1
            print(f'FAILED: {test_fn.__name__}  —  {exc}')
            import traceback
            traceback.print_exc()

    print('\n' + '='*70)
    print(f'ExoMol API test results: {passed} passed, {failed} failed out of {len(tests)}')
    print('='*70)
    sys.exit(1 if failed else 0)
