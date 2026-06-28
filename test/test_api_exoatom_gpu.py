"""
Test all PyExoCross API functions using ExoAtom Ar parameters.

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
    read_path='/public/home/zhangjingxin/LHD/Program/Databases/ExoAtom/', # '/Users/beryl/Academic/UCL/PhD/Data/database/ExoAtom/', 
    save_path='/public/home/zhangjingxin/LHD/Program/Data/pyexocross/gpu/', # '/Users/beryl/Academic/UCL/PhD/Data/pyexocross/gpu/',
    logs_path='/public/home/zhangjingxin/LHD/Program/Data/pyexocross/gpu/log/Ar_ExoAtom_gpu.log', # '/Users/beryl/Academic/UCL/PhD/Data/pyexocross/gpu/log/Ar_ExoAtom_gpu.log',
    cache='parquet', 
)


# Cores and chunks
COMPUTE_PARAMS = dict(
    ncputrans=4,                    # Number of CPU threads for each transition file (default: 4)
    ncpufiles=1,                    # Number of CPU files for transition calculation (default: 1)
    chunk_size=10000,               # Chunk size for transition calculation (default: 100000)
    device='GPU',                   # Device: 'CPU' or 'GPU' (default: 'CPU')
    gpu_backend='AUTO',             # GPU backend: 'AUTO', 'CUDA', 'PyTorch-CUDA', 'CuPy-CUDA', 'MPS' (default: 'AUTO')
    # gpu_batch_lines=8192,         # GPU line-batch size (default: 8192)
    # gpu_batch_grid=256,           # GPU grid-batch size (default: 256)
)


# Spectral range parameters
RANGE_PARAMS = dict(
    wn_wl='WN',                     # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
    wn_wl_unit='cm-1',              # Unit for wavenumber (default: cm⁻¹)
    min_range=0,                    # Minimum wavenumber in unit of cm⁻¹
    max_range=115400,               # Maximum wavenumber in unit of cm⁻¹
    unc_filter=None,                # Uncertainty filter (default: None)
)


# NLTE parameters (needed by stick_spectra and cross_sections)
NLTE_PARAMS = dict(
    nlte_method='P',                # Non-LTE: Population 
    nlte_path='/public/home/zhangjingxin/LHD/Program/Databases/ExoAtom/Ar/NIST/Ar_Ids.csv', # '/home/jingxin/LHD/Program/Databases/ExoAtom/Ar/NIST/Ar_Ids.csv',
)


data=px.load(
    **COMMON,
    **COMPUTE_PARAMS,
    **RANGE_PARAMS,
)


# ---------------------------------------------------------------------------
# Test functions
# ---------------------------------------------------------------------------
def test_cooling_functions():
    """Test cooling function calculation."""
    print('\n' + '='*70)
    print('TEST: px.cooling_functions()')
    print('='*70)
    px.cooling_functions(
        data=data,
        ntemp=1,                     # Number of temperature steps in unit of K (default: 1)
        tmax=600,                   # Maximum temperature in unit of K (default: 5000)
    )
    print('PASSED: cooling_functions()')


def test_stick_spectra():
    """Test stick spectra calculation."""
    print('\n' + '='*70)
    print('TEST: px.stick_spectra()')
    print('='*70)
    px.stick_spectra(
        data=data,
        # **NLTE_PARAMS,            # If Non-LTE is enabled, this parameter is required.
        temperatures=[1000,2000],   # Temperature in unit of K
        abs_emi='Ab',               # Absorption or emission (default: 'Absorption')
        threshold=None,             # Threshold filter (default: None)
        qns_filter={
            'configuration': [],
            'Multiple': [],
            'parity': [],
        },
        plot=True,                  # Whether to plot results (default: False)
        plot_method='log',          # Plot in linear (lin) or logarithm (log) (default: 'log')
        plot_wn_wl='Wn',            # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
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
        temperatures=[1000,2000],   # Temperature in unit of K
        pressures=[0.5,1.0],        # Pressure in bar (default: [1.0])
        abs_emi='Ab',               # Absorption or emission (default: 'Absorption')
        threshold=None,             # Threshold filter (default: None)
        bin_size=0.1,               # Bin size for wavenumber grid
        profile='SciPyVoigt',       # Line profile name (default: 'Gaussian')
        predissociation=False,      # Predissociation (default: False)
        cutoff=10.0,                # Cutoff distance in cm⁻¹ (default: None)
        broadeners=['Default'],     # Broadening species (default: ['Default'])
        ratios=[1.0],               # Broadening ratios (default: [1.0])
        alpha_hwhm=None,            # Constant Doppler HWHM (None, will calculate from broadening) or custom value (default: 3.0)
        gamma_hwhm=None,            # Constant Lorentzian HWHM (None, will calculate from broadening) or custom value (default: 0.5)
        plot=True,                  # Whether to plot results (default: False)
        plot_method='log',          # Plot in linear (lin) or logarithm (log)
        plot_wn_wl='WN',            # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='cm-1',           # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,          # Lower limit for y-axis (default: 1e-30 cm²/molecule)
    )
    print('PASSED: cross_sections()')


def test_stick_spectra_cross_section():
    """Test stick spectra and cross section calculation."""
    print('\n' + '='*70)
    print('TEST: px.stick_spectra_cross_section()')
    print('='*70)
    px.stick_spectra_cross_section(
        data=data,
        # **NLTE_PARAMS,            # If Non-LTE is enabled, this parameter is required.
        temperatures=[1000,2000],   # Temperature in unit of K
        pressures=[0.5,1.0],        # Pressure in bar (default: [1.0])
        abs_emi='Ab',               # Absorption or emission (default: 'Absorption')
        threshold=None,             # Threshold filter (default: None)
        bin_size=0.1,               # Bin size for wavenumber grid
        profile='SciPyVoigt',       # Line profile name (default: 'Gaussian')
        predissociation=False,      # Predissociation (default: False)
        cutoff=10.0,                # Cutoff distance in cm⁻¹ (default: None)
        broadeners=['Default'],     # Broadening species (default: ['Default'])
        ratios=[1.0],               # Broadening ratios (default: [1.0])
        alpha_hwhm=None,            # Constant Doppler HWHM (None, will calculate from broadening) or custom value (default: 3.0)
        gamma_hwhm=None,            # Constant Lorentzian HWHM (None, will calculate from broadening) or custom value (default: 0.5)
        plot=True,                  # Whether to plot results (default: False)
        plot_method='log',          # Plot in linear (lin) or logarithm (log)
        plot_wn_wl='WN',            # Wavenumber (wn in unit cm⁻¹) or wavelength (wl in unit[nm or um]) (default: 'WN')
        plot_unit='cm-1',           # Unit for plotting axis (default: cm⁻¹)
        limit_yaxis=1e-30,          # Lower limit for y-axis (default: 1e-30 cm²/molecule)
    )
    print('PASSED: stick_spectra_cross_section()')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print('PyExoCross API Test — ExoAtom Ar')
    print(f'pyexocross version: {px.__version__}')

    tests = [
        test_cooling_functions,
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
    print(f'ExoAtom API test results: {passed} passed, {failed} failed out of {len(tests)}')
    print('='*70)
    sys.exit(1 if failed else 0)
