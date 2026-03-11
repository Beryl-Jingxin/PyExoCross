"""
Test all PyExoCross API functions using ExoMolHR NO parameters.

Parameters are derived from: .input/AlCl_ExoMolHR.inp
"""
import os
import sys

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
    database='ExoMolHR',
    molecule='NO',
    isotopologue='14N-16O',
    species_id=81,
    read_path='/Users/beryl/Academic/UCL/PhD/Data/database/ExoMolHR/',
    save_path='/Users/beryl/Academic/UCL/PhD/Data/pyexocross/',
    logs_path='/Users/beryl/Academic/UCL/PhD/Data/pyexocross/log/test_api_exomolhr_no.log',
)

# Spectral range parameters
RANGE_PARAMS = dict(
    temperatures=[296, 1000],   # Temperature in unit of K
    wn_wl='WN',                 # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
    wn_wl_unit='cm-1',          # Unit for wavenumber (default: cm⁻¹)
    min_range=24,               # Minimum wavenumber in unit of cm⁻¹
    max_range=53452,            # Maximum wavenumber in unit of cm⁻¹
    abs_emi='Ab',               # Absorption or emission (default: 'Absorption')
    unc_filter=0.01,            # Uncertainty filter (default: None)
    threshold=1e-30,            # Threshold filter (default: None)
)

# NLTE parameters (needed by stick_spectra and cross_sections)
NLTE_PARAMS = dict(
    nlte_method='T',                   # Non-LTE: 2T NLTE (default: 'T')
    tvib_list=[1000, 2000],            # Vibrational temperatures in unit of K
    trot_list=[300],                   # Rotational temperatures in unit of K
    vib_label=['v', 'ElecState'],      # Vibrational quantum numbers
    rot_label=['J', 'e/f'],            # Rotational quantum numbers
)

# Cores and chunks
COMPUTE_PARAMS = dict(
    ncputrans=1,                # Number of CPU threads for each transition file (default: 4)
    ncpufiles=1,                # Number of CPU files for transition calculation (default: 1)
    chunk_size=10000,           # Chunk size for transition calculation (default: 100000)
)

# ---------------------------------------------------------------------------
# Test functions
# ---------------------------------------------------------------------------
def test_conversion():
    """Test ExoMolHR -> HITRAN conversion."""
    print('\n' + '='*70)
    print('TEST: px.conversion()  [ExoMolHR -> HITRAN]')
    print('='*70)
    px.conversion(
        **COMMON,
        **COMPUTE_PARAMS,
        conversion_format='HITRAN',
        conversion_min_freq=24,         # Minimum wavenumber in unit of cm⁻¹
        conversion_max_freq=53452,      # Maximum wavenumber in unit of cm⁻¹
        conversion_unc=0.01,            # Uncertainty filter (default: None)
        conversion_threshold=1e-30,     # Threshold filter (default: None)
        global_qn_label_list=['ElecState', 'v', 'Omega'],       # Quantum number label for global quantum numbers
        global_qn_format_list=['%9s', '%2d', '%4s'],            # Quantum number format for global quantum numbers
        local_qn_label_list=['J', 'e/f'],                       # Quantum number label for local quantum numbers
        local_qn_format_list=['%7.1f', '%1s'],                  # Quantum number format for local quantum numbers
    )
    print('PASSED: conversion()')
    
def test_stick_spectra():
    """Test stick spectra calculation."""
    print('\n' + '='*70)
    print('TEST: px.stick_spectra()')
    print('='*70)
    px.stick_spectra(
        **COMMON,
        # **NLTE_PARAMS,              # If Non-LTE is enabled, this parameter is required.
        **RANGE_PARAMS,
        **COMPUTE_PARAMS,
        plot=True,                    # Whether to plot results (default: False)
        plot_method='log',            # Plot in linear (lin) or logarithm (log) (default: 'log')
        plot_wn_wl='WN',              # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='nm',               # Unit for plotting axis (default: nm)
        limit_yaxis=1e-30,            # Lower limit for y-axis (default: 1e-30 cm/molecule)
    )
    print('PASSED: stick_spectra()')
    
def test_cross_sections():
    """Test cross section calculation."""
    print('\n' + '='*70)
    print('TEST: px.cross_sections()')
    print('='*70)
    px.cross_sections(
        **COMMON,
        # **NLTE_PARAMS,                # If Non-LTE is enabled, this parameter is required.
        **RANGE_PARAMS,
        **COMPUTE_PARAMS,
        pressures=[1.0],                # Pressure in unit bar (default: [1.0])
        bin_size=1,                     # Bin size for wavenumber grid 
        profile='SciPyVoigt',           # Line profile name (default: 'Gaussian')
        predissociation=False,          # Predissociation (default: False)
        cutoff=25.0,                    # Cutoff distance in cm⁻¹ (default: None)
        broadeners=['Default'],         # Broadening species (default: ['Default'])
        ratios=[1.0],                   # Broadening ratios (default: [1.0])
        alpha_hwhm=3.0,                 # Constant Doppler HWHM (None, will calculate from broadening) or custom value (default: 3.0)
        gamma_hwhm=None,                # Constant Lorentzian HWHM (None, will calculate from broadening) or custom value (default: 0.5)
        plot=True,                      # Whether to plot results (default: False)
        plot_method='log',              # Plot in linear (lin) or logarithm (log) (default: 'log')
        plot_wn_wl='WN',                # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='cm-1',               # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,              # Lower limit for y-axis (default: 1e-30 cm²/molecule)
    )
    print('PASSED: cross_sections()')

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print('PyExoCross API Test — ExoMolHR NO')
    print(f'pyexocross version: {px.__version__}')

    tests = [
        test_conversion,
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
    print(f'ExoMolHR API test results: {passed} passed, {failed} failed out of {len(tests)}')
    print('='*70)
    sys.exit(1 if failed else 0)