"""
Test all PyExoCross API functions using ExoMol MgH parameters.

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
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/Users/beryl/Academic/UCL/PhD/Data/database/ExoMol/', #'/home/jingxin/LHD/Program/Databases/ExoMol/',
    save_path='/Users/beryl/Academic/UCL/PhD/Data/pyexocross/', #'/home/jingxin/LHD/Program/Data/pyexocross/',
    logs_path='/Users/beryl/Academic/UCL/PhD/Data/pyexocross/log/MgH_ExoMol_parquet.log', #'/home/jingxin/LHD/Program/Data/pyexocross/log/test_api_exomol.log',
    cache='parquet',
)


# Cores and chunks
COMPUTE_PARAMS = dict(
    ncputrans=4,                # Number of CPU threads for each transition file (default: 4)
    ncpufiles=1,                # Number of CPU files for transition calculation (default: 1)
    chunk_size=500000,          # Chunk size for transition calculation (default: 100000)
    device='CPU',               # Device: 'CPU' or 'GPU' (default: 'CPU')
)


# Spectral range parameters
RANGE_PARAMS = dict(
    wn_wl='WN',                         # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
    wn_wl_unit='cm-1',                  # Unit for wavenumber (default: cm⁻¹)
    min_range=0,                        # Minimum wavenumber in unit of cm⁻¹
    max_range=30000,                    # Maximum wavenumber in unit of cm⁻¹
    unc_filter=0.5,                     # Uncertainty filter (default: None)
)


# NLTE parameters (needed by stick_spectra and cross_sections)
NLTE_PARAMS = dict(
    nlte_method='T',                    # Non-LTE: 2T NLTE (default: 'T')
    tvib_list=[3000],                   # Vibrational temperatures in unit of K
    trot_list=[296],                    # Rotational temperatures in unit of K
    vib_label=['v', 'ElecState'],       # Vibrational quantum numbers
    rot_label=['J', 'e/f'],             # Rotational quantum numbers
)


data=px.load(
    **COMMON,
    **COMPUTE_PARAMS,
    **RANGE_PARAMS,
)


# ---------------------------------------------------------------------------
# Test functions
# ---------------------------------------------------------------------------
def test_download():
    """Test download ExoMol files."""
    print('\n' + '='*70)
    print('TEST: px.download()')
    print('='*70)
    px.download(    
        file_path='/Users/beryl/Academic/UCL/PhD/Data/database/ExoMol/',    # Write that file_path or save_path are the same in download
        database='ExoMol',
        species_info={
            'MgH': {
                '24Mg-1H': {'wn_range': None},
                '25Mg-1H': {'wn_range': None},
            },
        },
        download=True,
    )
    print('PASSED: download()')
    
    
def test_conversion():
    """Test ExoMol -> HITRAN conversion."""
    print('\n' + '='*70)
    print('TEST: px.conversion()  [ExoMol -> HITRAN]')
    print('='*70)
    px.conversion(
        data=data,
        conversion_format='HITRAN',
        conversion_min_freq=0,          # Minimum wavenumber in unit of cm⁻¹
        conversion_max_freq=30000,      # Maximum wavenumber in unit of cm⁻¹
        conversion_unc=None,            # Uncertainty filter (default: None)
        conversion_threshold=None,      # Threshold filter (default: None)
        global_qn_label_list=['ElecState', 'v', 'Omega'],      # Quantum number label for global quantum numbers
        global_qn_format_list=['%9s', '%2d', '%4s'],           # Quantum number format for global quantum numbers
        local_qn_label_list=['J', 'e/f'],                      # Quantum number label for local quantum numbers
        local_qn_format_list=['%5.1f', '%2s'],                 # Quantum number format for local quantum numbers
    )
    print('PASSED: conversion()')


def test_partition_functions():
    """Test partition function calculation."""
    print('\n' + '='*70)
    print('TEST: px.partition_functions()')
    print('='*70)
    px.partition_functions(
        data=data,
        ntemp=1,                     # Number of temperature steps in unit of K (default: 1)
        tmax=5000,                   # Maximum temperature in unit of K (default: 5000)
    ) 
    print('PASSED: partition_functions()')


def test_specific_heats():
    """Test specific heat calculation."""
    print('\n' + '='*70)
    print('TEST: px.specific_heats()')
    print('='*70)
    px.specific_heats(
        data=data,
        ntemp=1,                     # Number of temperature steps in unit of K (default: 1)
        tmax=5000,                   # Maximum temperature in unit of K (default: 5000)
    )
    print('PASSED: specific_heats()')


def test_cooling_functions():
    """Test cooling function calculation."""
    print('\n' + '='*70)
    print('TEST: px.cooling_functions()')
    print('='*70)
    px.cooling_functions(
        data=data,
        ntemp=1,                     # Number of temperature steps in unit of K (default: 1)
        tmax=5000,                   # Maximum temperature in unit of K (default: 5000)
    )
    print('PASSED: cooling_functions()')


def test_lifetimes():
    """Test radiative lifetime calculation."""
    print('\n' + '='*70)
    print('TEST: px.lifetimes()')
    print('='*70)
    px.lifetimes(
        data=data,
        compress=False,               # Whether to compress the states file (default: False)
    )
    print('PASSED: lifetimes()')


def test_oscillator_strengths():
    """Test oscillator strength calculation."""
    print('\n' + '='*70)
    print('TEST: px.oscillator_strengths()')
    print('='*70)
    px.oscillator_strengths(
        data=data,
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
        data=data,
        # **NLTE_PARAMS,            # If Non-LTE is enabled, this parameter is required.
        temperatures=[296, 3000],   # Temperature in unit of K
        abs_emi='Ab',               # Absorption or emission (default: 'Absorption')
        threshold=None,             # Threshold filter (default: None)
        # qns_filter={
        #     '+/-': [],
        #     'e/f': [],
        #     'ElecState': [],
        #     'v': ['0,', '1,', '2,', '3,', '4,', ',0', ',1', ',2', ',3', ',4'],
        # },
        plot=True,                  # Whether to plot results (default: False)
        plot_method='log',          # Plot in linear (lin) or logarithm (log) (default: 'log')
        plot_wn_wl='WN',            # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WL')
        plot_unit='cm-1',           # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,          # Lower limit for y-axis (default: 1e-30 cm/molecule)
    )
    print('PASSED: stick_spectra()')


def test_cross_sections():
    """Test cross section calculation."""
    print('\n' + '='*70)
    print('TEST: px.cross_sections()')
    print('='*70)
    px.cross_sections(
        data=data,
        # **NLTE_PARAMS,            # If Non-LTE is enabled, this parameter is required.
        temperatures=[296, 3000],   # Temperature in unit of K
        pressures=[1.0],            # Pressure in unit bar (default: [1.0])
        abs_emi='Ab',               # Absorption or emission (default: 'Absorption')
        threshold=None,             # Threshold filter (default: None)
        bin_size=0.1,               # Bin size for wavenumber grid 
        profile='Gaussian',         # Line profile name (default: 'Gaussian')
        predissociation=False,      # Predissociation (default: False)
        cutoff=25.0,                # Cutoff distance in cm⁻¹ (default: None)
        broadeners=['Default'],     # Broadening species (default: ['Default'])
        ratios=[1.0],               # Broadening ratios (default: [1.0])
        alpha_hwhm=3,               # Constant Doppler HWHM (None, will calculate from broadening) or custom value (default: 3.0)
        gamma_hwhm=None,            # Constant Lorentzian HWHM (None, will calculate from broadening) or custom value (default: 0.5)
        plot=True,                  # Whether to plot results (default: False)
        plot_method='log',          # Plot in linear (lin) or logarithm (log) (default: 'log')
        plot_wn_wl='WN',            # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='cm-1',           # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,          # Lower limit for y-axis (default: 1e-30 cm²/molecule)
    )
    print('PASSED: cross_sections()')
    
    
def test_stick_spectra_cross_section():
    """Test stick spectra and cross section simultaneous calculation."""
    print('\n' + '='*70)
    print('TEST: px.stick_spectra_cross_section()')
    print('='*70)
    px.stick_spectra_cross_section(
        data=data,
        # **NLTE_PARAMS,            # If Non-LTE is enabled, this parameter is required.
        temperatures=[296, 3000],   # Temperature in unit of K
        pressures=[1.0],            # Pressure in unit bar (default: [1.0])
        wn_wl='WN',                 # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        wn_wl_unit='cm-1',          # Unit for wavenumber (default: cm⁻¹)
        min_range=0,                # Minimum wavenumber in unit of cm⁻¹
        max_range=13190,            # Maximum wavenumber in unit of cm⁻¹
        abs_emi='Ab',               # Absorption or emission (default: 'Absorption')
        threshold=None,             # Threshold filter (default: None)
        bin_size=0.1,               # Bin size for wavenumber grid 
        profile='SciPyVoigt',       # Line profile name (default: 'Gaussian')
        predissociation=False,      # Predissociation (default: False)
        cutoff=25.0,                # Cutoff distance in cm⁻¹ (default: None)
        broadeners=['Default'],     # Broadening species (default: ['Default'])
        ratios=[1.0],               # Broadening ratios (default: [1.0])
        alpha_hwhm=3,               # Constant Doppler HWHM (None, will calculate from broadening) or custom value (default: 3.0)
        gamma_hwhm=None,            # Constant Lorentzian HWHM (None, will calculate from broadening) or custom value (default: 0.5)
        plot=True,                  # Whether to plot results (default: False)
        plot_method='log',          # Plot in linear (lin) or logarithm (log) (default: 'log')
        plot_wn_wl='Wn',            # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='cm-1',           # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,          # Lower limit for y-axis (default: 1e-30 cm²/molecule)
    )
    print('PASSED: stick_spectra_cross_section()')



# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print('PyExoCross API Test — ExoMol MgH')
    print(f'pyexocross version: {px.__version__}')

    tests = [
        test_download,
        test_conversion,
        test_partition_functions,
        test_specific_heats,
        test_cooling_functions,
        test_lifetimes,
        test_oscillator_strengths,
        test_stick_spectra,
        test_cross_sections,
        test_stick_spectra_cross_section,
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
