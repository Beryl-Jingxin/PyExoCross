"""
High-level API functions for PyExoCross.

Usage::

    import pyexocross as px

    # Run all functions specified in an .inp file
    px.run('/path/to/config.inp')

    # Or call individual functions with keyword arguments
    px.cross_sections(
        database='ExoMol',
        molecule='MgH',
        isotopologue='24Mg-1H',
        dataset='XAB',
        read_path='/path/to/databases/',
        save_path='/path/to/output/',
        logs_path='/path/to/output/log/run.log',
        temperatures=[1000, 2000],
        pressures=[1.0],
        profile='SciPyVoigt',
    )

Available functions (snake_case, following PEP 8):

- ``px.run(inp_filepath)``         -- Run all functions from .inp file
- ``px.conversion(...)``           -- ExoMol-to-HITRAN or HITRAN-to-ExoMol
- ``px.partition_functions(...)``   -- Partition function calculation
- ``px.specific_heats(...)``       -- Specific heat calculation
- ``px.cooling_functions(...)``    -- Cooling function calculation
- ``px.lifetimes(...)``            -- Radiative lifetime calculation
- ``px.oscillator_strengths(...)`` -- Oscillator strength calculation
- ``px.stick_spectra(...)``        -- Stick spectra calculation
- ``px.cross_sections(...)``       -- Cross section calculation
"""
import os
from pyexocross.config import Config
from pyexocross.core import get_results


def _ensure_logging(inp_filepath=None, logs_path=None):
    """
    Set up logging to file if not already active.

    Mirrors the terminal behaviour: all ``print`` output is duplicated to a
    dated log file via ``TeeStream``.  If logging has already been
    initialised (e.g. by a prior API call or by ``run.py``), this is a
    no-op so that a single log file collects the complete session.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to an .inp file from which the log path is parsed.
    logs_path : str, optional
        Explicit log file path (used when calling with keyword args).
    """
    import pyexocross.base.log as _log_mod

    # If logging is already active, nothing to do
    if _log_mod._LOG_FILE_HANDLE is not None:
        return

    if inp_filepath is not None:
        # Parse the LogFilePath row from the .inp file
        log_path = _log_mod.parse_logging_info(inp_filepath)
    elif logs_path:
        # Ensure log directory exists
        log_dir = os.path.dirname(logs_path)
        if log_dir:
            from pyexocross.base.utils import ensure_dir
            ensure_dir(log_dir + '/')
        log_path = logs_path
    else:
        # No path given -- skip automatic logging
        return

    _log_mod.setup_logging(log_path)


def _remap_plot_kwargs(kwargs, plot_map):
    """Translate generic plot kwarg names to function-specific ones."""
    for generic, specific in plot_map.items():
        if generic in kwargs:
            kwargs[specific] = kwargs.pop(generic)


# ---------------------------------------------------------------------------
# Run (all functions from .inp file)
# ---------------------------------------------------------------------------
def run(inp_filepath):
    """
    Run all enabled functions from an .inp configuration file.

    This is the programmatic equivalent of running
    ``python3 run.py -p <inp_filepath>`` from the terminal.  All functions
    that are enabled in the .inp file (Conversion, PartitionFunctions,
    StickSpectra, CrossSections, etc.) will be executed, and the full log
    output will be saved to the path specified by ``LogFilePath`` in the
    .inp file.

    Parameters
    ----------
    inp_filepath : str
        Path to the .inp configuration file.

    Examples
    --------
    >>> import pyexocross as px
    >>> px.run('/path/to/MgH_ExoMol.inp')
    """
    _ensure_logging(inp_filepath=inp_filepath)
    config = Config(inp_filepath=inp_filepath)
    get_results(config)


# ---------------------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------------------
def conversion(inp_filepath=None, **kwargs):
    """
    Convert between ExoMol and HITRAN line-list formats.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        database : str
            'ExoMol' (format 1) or 'HITRAN'/'HITEMP' (format 2).
        molecule : str
        isotopologue : str
        dataset : str
        read_path : str
        save_path : str, optional
        logs_path : str, optional
            Log file path. Output is automatically saved when provided.
        species_id : int, optional
            Required for HITRAN/HITEMP databases.
        conversion_format : int
            1 for ExoMol -> HITRAN (default), 2 for HITRAN -> ExoMol.
        conversion_min_freq : float, optional
            Minimum frequency for conversion (default: 0).
        conversion_max_freq : float, optional
            Maximum frequency for conversion (default: 30000).
        conversion_unc : float or None, optional
            Uncertainty filter. ``None`` disables (default: ``None``).
        conversion_threshold : float or None, optional
            Intensity threshold. ``None`` disables (default: ``None``).
        global_qn_label_list : list of str, optional
            Required for ExoMol -> HITRAN conversion.
        global_qn_format_list : list of str, optional
        local_qn_label_list : list of str, optional
        local_qn_format_list : list of str, optional
        qnslabel_list : list of str, optional
            Quantum number labels.
        qnsformat_list : list of str, optional
            Quantum number formats.
        ncputrans : int, optional
            CPU cores for transitions (default: 4).
        ncpufiles : int, optional
            CPU cores for file I/O (default: 1).
        chunk_size : int, optional
            Chunk size for transitions (default: 100000).
    """
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = Config(inp_filepath=inp_filepath, conversion=1, **kwargs)
    get_results(config)


# Keep legacy aliases for backward compatibility
def convert_exomol_to_hitran(inp_filepath=None, **kwargs):
    """Convert ExoMol to HITRAN format. See :func:`conversion`."""
    kwargs.setdefault('conversion_format', 1)
    conversion(inp_filepath=inp_filepath, **kwargs)


def convert_hitran_to_exomol(inp_filepath=None, **kwargs):
    """Convert HITRAN to ExoMol format. See :func:`conversion`."""
    kwargs.setdefault('conversion_format', 2)
    conversion(inp_filepath=inp_filepath, **kwargs)


# ---------------------------------------------------------------------------
# Partition functions
# ---------------------------------------------------------------------------
def partition_functions(inp_filepath=None, **kwargs):
    """
    Calculate partition functions.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        database : str
        molecule : str (or atom for ExoAtom)
        isotopologue : str
        dataset : str
        read_path : str
        save_path : str, optional
        logs_path : str, optional
        species_id : int, optional
        ntemp : int, optional
            Temperature step interval (default: 1).
        tmax : int, optional
            Maximum temperature in K (default: 5000).
        ncputrans : int, optional
            CPU cores for transitions (default: 4).
        ncpufiles : int, optional
            CPU cores for file I/O (default: 1).
        chunk_size : int, optional
            Chunk size for transitions (default: 100000).
    """
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = Config(inp_filepath=inp_filepath, partition_functions=1, **kwargs)
    get_results(config)


# Legacy alias
partition_function = partition_functions


# ---------------------------------------------------------------------------
# Specific heats
# ---------------------------------------------------------------------------
def specific_heats(inp_filepath=None, **kwargs):
    """
    Calculate specific heat capacities.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        Same as :func:`partition_functions`.
    """
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = Config(inp_filepath=inp_filepath, specific_heats=1, **kwargs)
    get_results(config)


# Legacy alias
specific_heat = specific_heats


# ---------------------------------------------------------------------------
# Cooling functions
# ---------------------------------------------------------------------------
def cooling_functions(inp_filepath=None, **kwargs):
    """
    Calculate cooling functions.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        Same as :func:`partition_functions`.
    """
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = Config(inp_filepath=inp_filepath, cooling_functions=1, **kwargs)
    get_results(config)


# Legacy alias
cooling_function = cooling_functions


# ---------------------------------------------------------------------------
# Lifetimes
# ---------------------------------------------------------------------------
def lifetimes(inp_filepath=None, **kwargs):
    """
    Calculate radiative lifetimes.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        database : str
        molecule : str (or atom for ExoAtom)
        isotopologue : str
        dataset : str
        read_path : str
        save_path : str, optional
        logs_path : str, optional
        species_id : int, optional
        compress : bool, optional
            ``True`` to save as .bz2, ``False`` for uncompressed
            (default: ``False``).
        ncputrans : int, optional
            CPU cores for transitions (default: 4).
        ncpufiles : int, optional
            CPU cores for file I/O (default: 1).
        chunk_size : int, optional
            Chunk size for transitions (default: 100000).
    """
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = Config(inp_filepath=inp_filepath, lifetimes=1, **kwargs)
    get_results(config)


# Legacy alias
lifetime = lifetimes


# ---------------------------------------------------------------------------
# Oscillator strengths
# ---------------------------------------------------------------------------
def oscillator_strengths(inp_filepath=None, **kwargs):
    """
    Calculate oscillator strengths.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        database : str
        molecule : str (or atom for ExoAtom)
        isotopologue : str
        dataset : str
        read_path : str
        save_path : str, optional
        logs_path : str, optional
        species_id : int, optional
        gf_or_f : str, optional
            'gf' for weighted oscillator strength, 'f' for f-value
            (default: 'f').
        ncputrans : int, optional
            CPU cores for transitions (default: 4).
        ncpufiles : int, optional
            CPU cores for file I/O (default: 1).
        chunk_size : int, optional
            Chunk size for transitions (default: 100000).
        plot : bool, optional
            Whether to plot results (default: ``False``).
        plot_method : str, optional
            'log' or 'linear' (default: 'log').
        plot_wn_wl : str, optional
            'WN' for wavenumber, 'WL' for wavelength (default: 'WN').
        plot_unit : str, optional
            Unit for plotting axis (default: 'cm-1').
        limit_yaxis : float, optional
            Lower limit for y-axis (default: 1e-30).
    """
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    _remap_plot_kwargs(kwargs, {
        'plot': 'plot_oscillator_strength',
        'plot_method': 'plot_oscillator_strength_method',
        'plot_wn_wl': 'plot_oscillator_strength_wn_wl',
        'plot_unit': 'plot_oscillator_strength_unit',
        'limit_yaxis': 'limit_yaxis_os',
    })
    config = Config(inp_filepath=inp_filepath, oscillator_strengths=1, **kwargs)
    get_results(config)


# Legacy alias
oscillator_strength = oscillator_strengths


# ---------------------------------------------------------------------------
# Stick spectra
# ---------------------------------------------------------------------------
def stick_spectra(inp_filepath=None, **kwargs):
    """
    Calculate stick spectra.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        database : str
            'ExoMol', 'ExoAtom', 'HITRAN', or 'HITEMP'.
        molecule : str
        isotopologue : str
        dataset : str
        read_path : str
        save_path : str, optional
        logs_path : str, optional
        species_id : int, optional
            Required for HITRAN/HITEMP.
        temperatures : list of float, optional
            Temperature(s) in K (default: [1000]).
        nlte_method : str, optional
            'L' for LTE (default), 'T' for two-temperature NLTE,
            'D' for density NLTE, 'P' for population NLTE.
        tvib_list : list of float, optional
            Vibrational temperatures for NLTE 'T'.
        trot_list : list of float, optional
            Rotational temperatures for NLTE 'T' or 'D'.
        vib_label : list of str, optional
            Vibrational quantum number labels for NLTE 'T'.
        rot_label : list of str, optional
            Rotational quantum number labels for NLTE 'T'.
        nlte_path : str, optional
            Path to NLTE data file (required for methods 'D' and 'P').
        min_range : float, optional
            Minimum wavenumber/wavelength (default: 0).
        max_range : float, optional
            Maximum wavenumber/wavelength (default: 30000).
        wn_wl : str, optional
            'WN' for wavenumber, 'WL' for wavelength (default: 'WN').
        wn_wl_unit : str, optional
            'cm-1', 'um', or 'nm' (default: 'cm-1').
        abs_emi : str, optional
            'Ab' for absorption, 'Em' for emission (default: 'Ab').
        threshold : float or None, optional
            Intensity threshold. ``None`` disables (default: ``None``).
        unc_filter : float or None, optional
            Uncertainty filter. ``None`` disables (default: ``None``).
        qnslabel_list : list of str, optional
            Quantum number labels.
        qnsformat_list : list of str, optional
            Quantum number formats.
        qns_filter : dict, optional
            Quantum number filter. Keys are QN labels, values are lists of
            accepted values (or ``None`` to accept all).
            Example: ``{'v': ['0,', '1,', '2,'], 'par': None}``.
        ncputrans : int, optional
            CPU cores for transitions (default: 4).
        ncpufiles : int, optional
            CPU cores for file I/O (default: 1).
        chunk_size : int, optional
            Chunk size for transitions (default: 100000).
        plot : bool, optional
            Whether to plot results (default: ``False``).
        plot_method : str, optional
            'log' or 'linear' (default: 'log').
        plot_wn_wl : str, optional
            'WN' for wavenumber, 'WL' for wavelength (default: 'WN').
        plot_unit : str, optional
            Unit for plotting axis (default: 'cm-1').
        limit_yaxis : float, optional
            Lower limit for y-axis (default: 1e-30).

    Examples
    --------
    >>> import pyexocross as px
    >>> px.stick_spectra(
    ...     database='ExoMol',
    ...     molecule='MgH',
    ...     isotopologue='24Mg-1H',
    ...     dataset='XAB',
    ...     read_path='/path/to/databases/',
    ...     temperatures=[1000, 2000],
    ... )
    """
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    _remap_plot_kwargs(kwargs, {
        'plot': 'plot_stick_spectra',
        'plot_method': 'plot_stick_spectra_method',
        'plot_wn_wl': 'plot_stick_spectra_wn_wl',
        'plot_unit': 'plot_stick_spectra_unit',
        'limit_yaxis': 'limit_yaxis_stick_spectra',
    })
    config = Config(inp_filepath=inp_filepath, stick_spectra=1, **kwargs)
    get_results(config)


# ---------------------------------------------------------------------------
# Cross sections
# ---------------------------------------------------------------------------
def cross_sections(inp_filepath=None, **kwargs):
    """
    Calculate cross sections.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        All parameters from :func:`stick_spectra`, plus:

        pressures : list of float, optional
            Pressure(s) in bar (default: [1.0]).
        profile : str, optional
            Line profile name (default: 'Gaussian').
            Choices include 'Gaussian', 'Doppler', 'Lorentzian',
            'Voigt', 'SciPyVoigt', 'PseudoVoigt', etc.
        bin_size : float, optional
            Bin size for wavenumber grid (default: 0.1 cm-1).
            Mutually exclusive with ``n_point``.
        n_point : int, optional
            Number of grid points. Mutually exclusive with ``bin_size``.
        cutoff : float or None, optional
            Cutoff distance in cm-1. ``None`` disables (default: ``None``).
        broadeners : list of str, optional
            Broadening species (default: ['Default']).
        ratios : list of float, optional
            Broadening ratios (default: [1.0]).
        predissociation : bool, optional
            ``True`` to include predissociation (default: ``False``).
        alpha_hwhm : float or None, optional
            Constant Doppler HWHM value. ``None`` to calculate from
            broadening parameters (default: 3.0).
        gamma_hwhm : float or None, optional
            Constant Lorentzian HWHM value. ``None`` to calculate from
            broadening parameters (default: ``None``).
        plot : bool, optional
            Whether to plot results (default: ``False``).
        plot_method : str, optional
            'log' or 'linear' (default: 'log').
        plot_wn_wl : str, optional
            'WN' for wavenumber, 'WL' for wavelength (default: 'WN').
        plot_unit : str, optional
            Unit for plotting axis (default: 'cm-1').
        limit_yaxis : float, optional
            Lower limit for y-axis (default: 1e-30).

    Examples
    --------
    >>> import pyexocross as px
    >>> px.cross_sections(
    ...     database='ExoMol',
    ...     molecule='MgH',
    ...     isotopologue='24Mg-1H',
    ...     dataset='XAB',
    ...     read_path='/path/to/databases/',
    ...     temperatures=[1000],
    ...     pressures=[1.0],
    ...     profile='SciPyVoigt',
    ... )
    """
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    _remap_plot_kwargs(kwargs, {
        'plot': 'plot_cross_section',
        'plot_method': 'plot_cross_section_method',
        'plot_wn_wl': 'plot_cross_section_wn_wl',
        'plot_unit': 'plot_cross_section_unit',
        'limit_yaxis': 'limit_yaxis_xsec',
    })
    config = Config(inp_filepath=inp_filepath, cross_sections=1, **kwargs)
    get_results(config)


# Legacy alias
cross_section = cross_sections
