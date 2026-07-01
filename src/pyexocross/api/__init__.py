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

- ``px.run(inp_filepath)``                  -- Run all functions from .inp file
- ``px.load(...)``                          -- Load reusable ExoMol/ExoAtom data
- ``px.conversion(...)``                    -- ExoMol-to-HITRAN or HITRAN-to-ExoMol
- ``px.partition_functions(...)``           -- Partition function calculation
- ``px.specific_heats(...)``                -- Specific heat calculation
- ``px.cooling_functions(...)``             -- Cooling function calculation
- ``px.lifetimes(...)``                     -- Radiative lifetime calculation
- ``px.oscillator_strengths(...)``          -- Oscillator strength calculation
- ``px.stick_spectra(...)``                 -- Stick spectra calculation
- ``px.cross_sections(...)``                -- Cross section calculation
- ``px.stick_spectra_cross_section(...)``   -- Stick spectra and Cross section calculation simultaneously
"""
import os
from ..config import Config
from ..core import get_results, printdatabaseinfo, printdeviceinfo
from ..base.log import close_logging
from ..base.utils import Timer
from ..database.data import LoadedData, loaddata


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


def _calculation_config(inp_filepath, data, kwargs, **flags):
    """Build a direct or LoadedData-backed calculation configuration."""
    if data is None:
        options = dict(kwargs)
        options.update(flags)
        return Config(inp_filepath=inp_filepath, **options)
    if inp_filepath is not None:
        raise ValueError('inp_filepath cannot be combined with data from px.load.')
    if not isinstance(data, LoadedData):
        raise TypeError('data must be returned by px.load.')
    config = data.configfor(**kwargs)
    for name, value in flags.items():
        setattr(config, name, value)
    validateloadrange(data, config)
    return config


def validateloadrange(data, config):
    """Reject calculation ranges that exceed the interval supplied to px.load."""
    if config.conversion == 1:
        requestedmin = config.conversion_min_freq
        requestedmax = config.conversion_max_freq
    elif config.stick_spectra == 1 or config.cross_sections == 1:
        requestedmin = config.min_wn
        requestedmax = config.max_wn
    else:
        return

    loadedmin = data.config.min_wn
    loadedmax = data.config.max_wn
    requiredmin = float(requestedmin)
    requiredmax = float(requestedmax)
    if requiredmin < loadedmin or requiredmax > loadedmax:
        recommendedmin = min(float(loadedmin), requiredmin)
        recommendedmax = max(float(loadedmax), requiredmax)
        if data.config.wn_wl == 'WL':
            unitfactor = 1e4 if data.config.wn_wl_unit == 'um' else 1e7
            loadmin = unitfactor / recommendedmax
            loadmax = unitfactor / recommendedmin if recommendedmin > 0 else float('inf')
            recommendation = (
                f'min_range={loadmin:g}, max_range={loadmax:g} '
                f'({data.config.wn_wl_unit})'
            )
        else:
            recommendation = (
                f'min_range={recommendedmin:g}, max_range={recommendedmax:g}'
            )
        raise ValueError(
            f'Requested calculation range {float(requestedmin):g}-{float(requestedmax):g} '
            f'cm-1 is not covered by px.load range '
            f'{float(loadedmin):g}-{float(loadedmax):g} cm-1. '
            f'Call px.load({recommendation}) again.'
        )


def requirealltransitions(data):
    """Expand range-loaded ExoMol/ExoAtom data when a whole-list calculation needs it."""
    if (
        data is None
        or data.config.database not in ('ExoMol', 'ExoAtom')
        or data.alltrans
    ):
        return data

    print('Loading all transitions required by this calculation ...')
    expanded = loaddata(
        data.config,
        cache=data.cache,
        cachedir=data.cachedir,
        maxmemory=data.config.max_memory,
        refresh=data.config.refresh_cache,
        alltrans=True,
        preparestates=data.preparedstates is not None,
    )
    data.states = expanded.states
    data.preparedstates = expanded.preparedstates
    data.transitions = expanded.transitions
    data.storage = expanded.storage
    data.alltrans = True
    return data


# ---------------------------------------------------------------------------
# Reusable data loading
# ---------------------------------------------------------------------------
def load(
    inp_filepath=None,
    cache='auto',
    cache_dir=None,
    max_memory=512,
    refresh=False,
    all_transitions=False,
    **kwargs,
):
    """
    Read and preprocess reusable line-list data.

    Parameters
    ----------
    cache : {'auto', 'parquet', 'none'}
        ``auto`` keeps small transitions in memory and converts large inputs to
        Parquet. ``parquet`` always uses a persistent Parquet cache. ``none``
        keeps the original transition files as streaming sources.
    cache_dir : str, optional
        Parquet cache directory. Defaults to ``.pyexocross_cache`` inside the
        selected input dataset directory.
    max_memory : int, optional
        Maximum estimated transition data retained by ``auto`` in memory, in
        MB. Default is 512.
    refresh : bool, optional
        Rebuild matching Parquet cache files.
    all_transitions : bool, optional
        Eagerly load every transition file. ExoMol/ExoAtom data is expanded
        automatically when ``lifetimes``, ``cooling_functions``, or
        ``oscillator_strengths`` first needs it. Default is False.
    **kwargs
        The same database, range, filtering, and preprocessing options accepted
        by :func:`stick_spectra`. ExoMol, ExoAtom, ExoMolHR, HITRAN, and HITEMP
        are supported.
    """
    if 'abundance' in kwargs:
        raise ValueError(
            'abundance is a calculation parameter. Pass it to stick_spectra, '
            'cross_section, or stick_spectra_cross_section instead of px.load.'
        )
    if 'cutoff' in kwargs:
        raise ValueError(
            'cutoff is a calculation parameter. Pass it to cross_section or '
            'stick_spectra_cross_section instead of px.load.'
        )
    suppliedplotkeys = sorted(
        key
        for key in kwargs
        if key == 'plot' or key.startswith('plot_') or key.startswith('limit_yaxis')
    )
    if suppliedplotkeys:
        raise ValueError(
            ', '.join(suppliedplotkeys)
            + ' are calculation parameters. Pass them to the corresponding '
            'calculation function instead of px.load.'
        )
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = Config(
        inp_filepath=inp_filepath,
        cache=cache,
        cache_dir=cache_dir,
        max_memory=max_memory,
        refresh=refresh,
        **kwargs,
    )
    config.to_globals()
    printdeviceinfo(config)
    printdatabaseinfo(config)
    print()
    timer = Timer().start()
    print('Loading reusable line-list data ...')
    data = loaddata(
        config,
        cache=cache,
        cachedir=cache_dir,
        maxmemory=max_memory,
        refresh=refresh,
        alltrans=all_transitions,
    )
    print('Finished loading reusable line-list data!')
    timer.end()
    return data


load_data = load


# ---------------------------------------------------------------------------
# Run (all functions from .inp file)
# ---------------------------------------------------------------------------
def run(inp_filepath, force_reload=False):
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
    force_reload : bool, optional
        If True, force re-parse of ``inp_filepath`` even when a cached
        configuration exists in this Python process. Default is False.

    Examples
    --------
    >>> import pyexocross as px
    >>> px.run('/path/to/MgH_ExoMol.inp')
    """
    _ensure_logging(inp_filepath=inp_filepath)
    config = Config(inp_filepath=inp_filepath, force_reload=force_reload)
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
            Maximum frequency for conversion (default: 1e10).
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
    data = kwargs.pop('data', None)
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = _calculation_config(inp_filepath, data, kwargs, conversion=1)
    get_results(config, data=data)


# Keep legacy aliases for backward compatibility
def convert_exomol_to_hitran(inp_filepath=None, **kwargs):
    """Convert ExoMol to HITRAN format. See :func:`conversion`."""
    kwargs.setdefault('conversion_format', 'HITRAN')
    conversion(inp_filepath=inp_filepath, **kwargs)
    
def convert_exomolhr_to_hitran(inp_filepath=None, **kwargs):
    """Convert ExoMolHR to HITRAN format. See :func:`conversion`."""
    kwargs.setdefault('conversion_format', 'HITRAN')
    conversion(inp_filepath=inp_filepath, **kwargs)
    
def convert_exoatom_to_hitran(inp_filepath=None, **kwargs):
    """Convert ExoAtom to HITRAN format. See :func:`conversion`."""
    kwargs.setdefault('conversion_format', 'HITRAN')
    conversion(inp_filepath=inp_filepath, **kwargs)

def convert_hitran_to_exomol(inp_filepath=None, **kwargs):
    """Convert HITRAN to ExoMol format. See :func:`conversion`."""
    kwargs.setdefault('conversion_format', 'ExoMol')
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
    data = kwargs.pop('data', None)
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = _calculation_config(
        inp_filepath,
        data,
        kwargs,
        partition_functions=1,
    )
    get_results(config, data=data)


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
    data = kwargs.pop('data', None)
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = _calculation_config(
        inp_filepath,
        data,
        kwargs,
        specific_heats=1,
    )
    get_results(config, data=data)


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
    data = kwargs.pop('data', None)
    data = requirealltransitions(data)
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = _calculation_config(inp_filepath, data, kwargs, cooling_functions=1)
    get_results(config, data=data)


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
    data = kwargs.pop('data', None)
    data = requirealltransitions(data)
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    config = _calculation_config(inp_filepath, data, kwargs, lifetimes=1)
    get_results(config, data=data)


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
    data = kwargs.pop('data', None)
    data = requirealltransitions(data)
    config = _calculation_config(inp_filepath, data, kwargs, oscillator_strengths=1)
    get_results(config, data=data)


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
        abundance : float, optional
            Isotopic abundance multiplier (default: 1.0).
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

        Notes
        -----
        ``wn_wl`` and ``wn_wl_unit`` control the calculation range, output file
        naming, and the first column in the saved stick spectra. Plotting
        options only control the figure x-axis.

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
    data = kwargs.pop('data', None)
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    _remap_plot_kwargs(kwargs, {
        'plot': 'plot_stick_spectra',
        'plot_method': 'plot_stick_spectra_method',
        'plot_wn_wl': 'plot_stick_spectra_wn_wl',
        'plot_unit': 'plot_stick_spectra_unit',
        'limit_yaxis': 'limit_yaxis_stick_spectra',
    })
    config = _calculation_config(
        inp_filepath,
        data,
        kwargs,
        stick_spectra=1,
    )
    get_results(config, data=data)


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
            Bin size in the selected ``wn_wl_unit`` (default: 0.1).
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
            Constant Doppler HWHM value. ``None`` calculates it from the
            molecular mass and temperature (default: ``None``).
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

        Notes
        -----
        ``wn_wl`` and ``wn_wl_unit`` control the calculation range, output file
        naming, cross-section grid unit, and the first column in the saved
        cross sections. With ``wn_wl='WL'``, ``bin_size`` is a wavelength
        interval in ``wn_wl_unit``. Line-profile widths and ``cutoff`` remain
        in cm-1 because profiles are evaluated internally in wavenumber space.

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
    data = kwargs.pop('data', None)
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    _remap_plot_kwargs(kwargs, {
        'plot': 'plot_cross_section',
        'plot_method': 'plot_cross_section_method',
        'plot_wn_wl': 'plot_cross_section_wn_wl',
        'plot_unit': 'plot_cross_section_unit',
        'limit_yaxis': 'limit_yaxis_xsec',
    })
    config = _calculation_config(
        inp_filepath,
        data,
        kwargs,
        cross_sections=1,
    )
    get_results(config, data=data)


# Legacy alias
cross_section = cross_sections


# ---------------------------------------------------------------------------
# Stick spectra and Cross sections simultaneously
# ---------------------------------------------------------------------------
def stick_spectra_cross_section(inp_filepath=None, **kwargs):
    """
    Calculate stick spectra and cross sections simultaneously.

    Parameters
    ----------
    inp_filepath : str, optional
        Path to .inp configuration file.
    **kwargs
        All parameters from :func:`stick_spectra` and :func:`cross_sections`.

        Notes
        -----
        The same ``wn_wl`` and ``wn_wl_unit`` selection is used for both saved
        outputs. The first columns of ``.stick`` and ``.xsec`` are wavenumber
        for ``wn_wl='WN'`` and wavelength for ``wn_wl='WL'``. Cross-section
        ``bin_size`` follows the same selected unit.

    Examples
    --------
    >>> import pyexocross as px
    >>> px.stick_spectra_cross_section(
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
    data = kwargs.pop('data', None)
    _ensure_logging(inp_filepath, kwargs.get('logs_path'))
    
    if 'plot' in kwargs:
        plot_val = kwargs.pop('plot')
        kwargs.setdefault('plot_stick_spectra', plot_val)
        kwargs.setdefault('plot_cross_section', plot_val)
    if 'plot_method' in kwargs:
        method_val = kwargs.pop('plot_method')
        kwargs.setdefault('plot_stick_spectra_method', method_val)
        kwargs.setdefault('plot_cross_section_method', method_val)
    if 'plot_wn_wl' in kwargs:
        wn_wl_val = kwargs.pop('plot_wn_wl')
        kwargs.setdefault('plot_stick_spectra_wn_wl', wn_wl_val)
        kwargs.setdefault('plot_cross_section_wn_wl', wn_wl_val)
    if 'plot_unit' in kwargs:
        unit_val = kwargs.pop('plot_unit')
        kwargs.setdefault('plot_stick_spectra_unit', unit_val)
        kwargs.setdefault('plot_cross_section_unit', unit_val)
    if 'limit_yaxis' in kwargs:
        limit_val = kwargs.pop('limit_yaxis')
        kwargs.setdefault('limit_yaxis_stick_spectra', limit_val)
        kwargs.setdefault('limit_yaxis_xsec', limit_val)

    config = _calculation_config(
        inp_filepath,
        data,
        kwargs,
        stick_spectra=1,
        cross_sections=1,
    )
    get_results(config, data=data)
