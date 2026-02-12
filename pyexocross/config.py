"""
Configuration management for PyExoCross.

Supports two input methods:
1. Reading from .inp configuration file
2. Direct parameter passing via function arguments
"""
import os
import re
import sys
import warnings
import numpy as np
import pandas as pd
from tabulate import tabulate

# Add src directory to path for imports
_current_dir = os.path.dirname(os.path.abspath(__file__))
# _current_dir is pyexocross/, so parent is project root
_project_root = os.path.dirname(_current_dir)
_src_dir = os.path.join(_project_root, 'src')
if _src_dir not in sys.path and os.path.exists(_src_dir):
    sys.path.insert(0, _src_dir)

from base.input import inp_para, parse_TP_values
from base.utils import ensure_dir
from base.constants import *

class Config:
    """
    Configuration class that manages all PyExoCross parameters.
    
    Can be initialized from an .inp file or from keyword arguments.
    """
    
    def __init__(self, inp_filepath=None, **kwargs):
        """
        Initialize configuration from .inp file or keyword arguments.

        Parameters
        ----------
        inp_filepath : str, optional
            Path to .inp configuration file. If provided, parameters are read from file.
        **kwargs
            Direct parameter values. Used when inp_filepath is None or to override file values.
        """
        if inp_filepath is not None:
            self._load_from_file(inp_filepath)
        
        # Always call _load_from_kwargs to handle overrides and defaults.
        # If inp_filepath was provided, kwargs will override file values.
        # If inp_filepath was None, kwargs (plus defaults) will set initial values.
        if inp_filepath is None or kwargs:
            self._load_from_kwargs(**kwargs)
        
        # Store on ConfigManager so get_config() works for both pathways.
        # Import from both possible module paths to handle Python's dual-path
        # module registration (base.config_manager vs src.base.config_manager).
        from base.config_manager import ConfigManager
        ConfigManager._last_config = self
        try:
            from src.base.config_manager import ConfigManager as _CM2
            _CM2._last_config = self
        except ImportError:
            pass
    
    def _load_from_file(self, inp_filepath):
        """Load configuration from .inp file."""
        # Use ConfigManager to cache parsed configurations
        from base.config_manager import ConfigManager
        params = ConfigManager.get_config(inp_filepath)
        self._set_attributes(params)
    
    def _load_from_kwargs(self, **kwargs):
        """Load configuration from keyword arguments with defaults."""
        # Database info
        self.database = kwargs.get('database', getattr(self, 'database', 'ExoMol'))
        self.molecule = kwargs.get('molecule', getattr(self, 'molecule', None))
        self.isotopologue = kwargs.get('isotopologue', getattr(self, 'isotopologue', None))
        self.dataset = kwargs.get('dataset', getattr(self, 'dataset', None))
        self.atom = kwargs.get('atom', getattr(self, 'atom', None))  # For ExoAtom
        self.species_id = kwargs.get('species_id', getattr(self, 'species_id', None))
        
        # File paths
        self.read_path = kwargs.get('read_path', getattr(self, 'read_path', './'))
        self.save_path = kwargs.get('save_path', getattr(self, 'save_path', './output/'))
        self.logs_path = kwargs.get('logs_path', getattr(self, 'logs_path', './pyexocross.log'))
        
        # Ensure directories exist
        ensure_dir(self.save_path)
        log_dir = os.path.dirname(self.logs_path)
        if log_dir:
            ensure_dir(log_dir + '/')
        
        # Function flags
        self.conversion = kwargs.get('conversion', getattr(self, 'conversion', 0))
        self.partition_functions = kwargs.get('partition_functions', getattr(self, 'partition_functions', 0))
        self.specific_heats = kwargs.get('specific_heats', getattr(self, 'specific_heats', 0))
        self.lifetimes = kwargs.get('lifetimes', getattr(self, 'lifetimes', 0))
        self.cooling_functions = kwargs.get('cooling_functions', getattr(self, 'cooling_functions', 0))
        self.oscillator_strengths = kwargs.get('oscillator_strengths', getattr(self, 'oscillator_strengths', 0))
        self.stick_spectra = kwargs.get('stick_spectra', getattr(self, 'stick_spectra', 0))
        self.cross_sections = kwargs.get('cross_sections', getattr(self, 'cross_sections', 0))
        # For ExoMol / ExoAtom, SpeciesID is only required for conversion
        # because HITRAN output needs molecule/isotopologue IDs.
        if (
            self.database in ('ExoMol', 'ExoAtom')
            and self.conversion == 1
            and (self.species_id is None or int(self.species_id) == 0)
        ):
            raise ValueError(
                "species_id is required when conversion is enabled for ExoMol/ExoAtom."
            )
        
        # Computational parameters
        self.ncputrans = kwargs.get('ncputrans', getattr(self, 'ncputrans', 4))
        self.ncpufiles = kwargs.get('ncpufiles', getattr(self, 'ncpufiles', 1))
        self.chunk_size = kwargs.get('chunk_size', getattr(self, 'chunk_size', 100000))
        
        # Temperature and pressure
        if 'temperatures' in kwargs:
            T_input = kwargs['temperatures']
            self.T_list = parse_TP_values(str(T_input[0])) if isinstance(T_input, (list, np.ndarray)) and len(T_input) == 1 else T_input
        elif not hasattr(self, 'T_list'):
            self.T_list = [1000]

        if 'pressures' in kwargs:
            P_input = kwargs['pressures']
            self.P_list = parse_TP_values(str(P_input[0])) if isinstance(P_input, (list, np.ndarray)) and len(P_input) == 1 else P_input
        elif not hasattr(self, 'P_list'):
            self.P_list = [1.0]
        
        # Wavenumber/wavelength settings
        self.wn_wl = kwargs.get('wn_wl', getattr(self, 'wn_wl', 'WN'))
        self.wn_wl_unit = kwargs.get('wn_wl_unit', getattr(self, 'wn_wl_unit', 'cm-1'))
        
        min_range = kwargs.get('min_range', getattr(self, 'min_wnl', 0))
        max_range = kwargs.get('max_range', getattr(self, 'max_wnl', 30000))
        
        self.min_wn = min_range if self.wn_wl == 'WN' else 1e4/max_range if self.wn_wl_unit == 'um' else 1e7/max_range
        self.max_wn = max_range if self.wn_wl == 'WN' else 1e4/min_range if self.wn_wl_unit == 'um' else 1e7/min_range
        self.min_wnl = min_range
        self.max_wnl = max_range
        
        # Grid settings (n_point and bin_size are mutually exclusive)
        _n_point = kwargs.get('n_point', None)
        _bin_size = kwargs.get('bin_size', None)
        if _n_point is not None and _bin_size is not None:
            raise ValueError("Provide either n_point or bin_size, not both.")
        
        if _n_point is not None:
            self.N_point = int(_n_point)
            self.bin_size = (self.max_wn - self.min_wn) / max(self.N_point - 1, 1)
        elif _bin_size is not None:
            bs = _bin_size if self.wn_wl == 'WN' else (_bin_size * 1e-4 if self.wn_wl_unit == 'um' else _bin_size * 1e-7)
            self.bin_size = bs
            self.N_point = int((self.max_wn - self.min_wn) / self.bin_size) + 1
        elif not hasattr(self, 'N_point') or not hasattr(self, 'bin_size'):
            # Default: bin_size = 0.1 cm-1
            self.bin_size = 0.1
            self.N_point = int((self.max_wn - self.min_wn) / self.bin_size) + 1
        self.wn_grid = np.linspace(self.min_wn, self.max_wn, self.N_point)
        
        # Filters (None = disabled, float/int value = enabled)
        current_threshold = getattr(self, 'threshold', 'None')
        _threshold = kwargs.get('threshold', current_threshold if current_threshold != 'None' else None)
        self.threshold = 'None' if _threshold is None else _threshold
        
        current_unc = getattr(self, 'unc_filter', 'None')
        _unc_filter = kwargs.get('unc_filter', current_unc if current_unc != 'None' else None)
        self.unc_filter = 'None' if _unc_filter is None else _unc_filter
        
        # NLTE settings
        self.nlte_method = kwargs.get('nlte_method', getattr(self, 'nlte_method', 'L'))
        self.tvib_list = list(kwargs.get('tvib_list', getattr(self, 'tvib_list', [])))
        self.trot_list = list(kwargs.get('trot_list', getattr(self, 'trot_list', [])))
        self.vib_label = kwargs.get('vib_label', getattr(self, 'vib_label', []))
        self.rot_label = kwargs.get('rot_label', getattr(self, 'rot_label', []))
        self.nlte_path = kwargs.get('nlte_path', getattr(self, 'nlte_path', None))
        
        # Line profile (uppercase and strip 'PRO' suffix to match inp_para convention)
        _profile_raw = kwargs.get('profile', getattr(self, 'profile', 'Gaussian'))
        self.profile = _profile_raw.upper().replace('PRO', '')
        
        _abs_emi_raw = kwargs.get('abs_emi', getattr(self, 'abs_emi', 'Ab'))
        self.abs_emi = 'Ab' if _abs_emi_raw.upper().startswith('AB') else 'Em'
        
        current_cutoff = getattr(self, 'cutoff', 'None')
        _cutoff = kwargs.get('cutoff', current_cutoff if current_cutoff != 'None' else None)
        self.cutoff = 'None' if _cutoff is None else _cutoff
        
        # Predissociation: bool (new) or legacy 'Y'/'N' string
        _predissoc_default = getattr(self, 'predissoc_yn', False)
        # Convert file-loaded 'Y'/'N' to bool for default logic
        if isinstance(_predissoc_default, str):
            _predissoc_default = (_predissoc_default == 'Y')
            
        _predissoc = kwargs.get('predissociation', kwargs.get('predissoc_yn', _predissoc_default))
        if isinstance(_predissoc, bool):
            self.predissoc_yn = 'Y' if _predissoc else 'N'
        else:
            self.predissoc_yn = _predissoc
        
        # Broadening
        self.broadeners = kwargs.get('broadeners', getattr(self, 'broadeners', ['Default']))
        self.ratios = kwargs.get('ratios', getattr(self, 'ratios', [1.0]))
        
        # Set data_info based on database
        if self.database == 'ExoAtom':
            if self.atom is None:
                raise ValueError("atom must be specified for ExoAtom database")
            self.data_info = [self.atom, self.dataset]
        else:
            if self.molecule is None or self.isotopologue is None or self.dataset is None:
                raise ValueError("molecule, isotopologue, and dataset must be specified")
            self.data_info = [self.molecule, self.isotopologue, self.dataset]
        
        # Set defaults for other parameters
        self._set_defaults(**kwargs)
        
        # Resolve database-specific metadata from definition files
        self._resolve_metadata()
        
        # NLTE mesh grid alignment (matching inp_para logic)
        if self.nlte_method == 'T':
            # Two-temperature NLTE: align Tvib_list x Trot_list
            if len(self.tvib_list) == 1 and len(self.trot_list) > 1:
                self.tvib_list = self.tvib_list * len(self.trot_list)
            elif len(self.tvib_list) > 1 and len(self.trot_list) == 1:
                self.trot_list = self.trot_list * len(self.tvib_list)
            elif len(self.tvib_list) > 1 and len(self.trot_list) > 1:
                grid_v, grid_r = np.meshgrid(self.tvib_list, self.trot_list, indexing='ij')
                self.tvib_list = grid_v.ravel().tolist()
                self.trot_list = grid_r.ravel().tolist()
            self.lte_nlte = '__2T-nlte'
        elif self.nlte_method == 'D':
            # Density NLTE: align T_list x Trot_list; requires nlte_path
            if not self.nlte_path:
                raise ValueError("nlte_path is required for NLTE method 'D'")
            self.vib_label = []
            self.rot_label = []
            self.tvib_list = []
            if len(self.T_list) == 1 and len(self.trot_list) > 1:
                self.T_list = self.T_list * len(self.trot_list)
            elif len(self.T_list) > 1 and len(self.trot_list) == 1:
                self.trot_list = self.trot_list * len(self.T_list)
            elif len(self.T_list) > 1 and len(self.trot_list) > 1:
                grid_t, grid_r = np.meshgrid(self.T_list, self.trot_list, indexing='ij')
                self.T_list = grid_t.ravel().tolist()
                self.trot_list = grid_r.ravel().tolist()
            self.lte_nlte = '__nvib-nlte'
        elif self.nlte_method == 'P':
            # Population NLTE: requires nlte_path; no extra temperature lists
            if not self.nlte_path:
                raise ValueError("nlte_path is required for NLTE method 'P'")
            self.vib_label = []
            self.rot_label = []
            self.tvib_list = []
            self.trot_list = []
            self.lte_nlte = '__pop-nlte'
        else:
            self.lte_nlte = ''
        
        # Derive photo suffix from predissoc_yn
        self.photo = '__photo' if self.predissoc_yn == 'Y' else ''
    
    def _set_defaults(self, **kwargs):
        """Set default values for parameters not specified, with kwargs overrides."""
        # Conversion defaults
        self.conversion_format = kwargs.get('conversion_format', getattr(self, 'conversion_format', 1))
        self.conversion_min_freq = kwargs.get('conversion_min_freq', getattr(self, 'conversion_min_freq', 0))
        self.conversion_max_freq = kwargs.get('conversion_max_freq', getattr(self, 'conversion_max_freq', 30000))
        
        curr_conv_unc = getattr(self, 'conversion_unc', 'None')
        _conv_unc = kwargs.get('conversion_unc', curr_conv_unc if curr_conv_unc != 'None' else None)
        self.conversion_unc = 'None' if _conv_unc is None else _conv_unc
        
        curr_conv_thr = getattr(self, 'conversion_threshold', 'None')
        _conv_thr = kwargs.get('conversion_threshold', curr_conv_thr if curr_conv_thr != 'None' else None)
        self.conversion_threshold = 'None' if _conv_thr is None else _conv_thr
        
        self.global_qn_label_list = kwargs.get('global_qn_label_list', getattr(self, 'global_qn_label_list', []))
        self.global_qn_format_list = kwargs.get('global_qn_format_list', getattr(self, 'global_qn_format_list', []))
        self.local_qn_label_list = kwargs.get('local_qn_label_list', getattr(self, 'local_qn_label_list', []))
        self.local_qn_format_list = kwargs.get('local_qn_format_list', getattr(self, 'local_qn_format_list', []))
        
        # Partition function defaults
        self.ntemp = kwargs.get('ntemp', getattr(self, 'ntemp', 1))
        self.tmax = kwargs.get('tmax', getattr(self, 'tmax', 5000))
        
        # Compress: bool (new) or legacy 'Y'/'N' string
        _compress_default = getattr(self, 'compress_yn', False)
        if isinstance(_compress_default, str):
             _compress_default = (_compress_default == 'Y')
             
        _compress = kwargs.get('compress', kwargs.get('compress_yn', _compress_default))
        if isinstance(_compress, bool):
            self.compress_yn = 'Y' if _compress else 'N'
        else:
            self.compress_yn = _compress

        self.gf_or_f = kwargs.get('gf_or_f', getattr(self, 'gf_or_f', 'f'))
        self.lte_nlte = ''
        self.photo = ''
        self.qnslabel_list = kwargs.get('qnslabel_list', getattr(self, 'qnslabel_list', []))
        self.qnsformat_list = kwargs.get('qnsformat_list', getattr(self, 'qnsformat_list', []))
        
        # QN Filter logic
        _qns_filter_raw = kwargs.get('qns_filter', None)
        
        if _qns_filter_raw is not None:
            # Explicit kwarg provided
            if isinstance(_qns_filter_raw, dict) and _qns_filter_raw:
                # Dict-based
                self.qns_label = list(_qns_filter_raw.keys())
                self.qns_value = [v if v else [''] for v in _qns_filter_raw.values()]
                self.qns_format = [
                    self.qnsformat_list[self.qnslabel_list.index(label)]
                    for label in self.qns_label
                    if label in self.qnslabel_list
                ]
                self.qns_filter = [
                    f"{label}[{';'.join(vals)}]" if vals != [''] else label
                    for label, vals in zip(self.qns_label, self.qns_value)
                ]
            else:
                # Legacy list provided in kwargs
                self.qns_label = kwargs.get('qns_label', [])
                self.qns_value = kwargs.get('qns_value', [])
                self.qns_format = kwargs.get('qns_format', [])
                self.qns_filter = _qns_filter_raw
        elif not hasattr(self, 'qns_filter'):
            # Not in kwargs, not in self -> default empty
            self.qns_label = []
            self.qns_value = []
            self.qns_format = []
            self.qns_filter = []
        # else: use existing self.qns_filter from file

        # Doppler HWHM
        _alpha = kwargs.get('alpha_hwhm', getattr(self, 'alpha_hwhm', 3.0))
        # Prioritize doppler_hwhm_yn if explicitly passed, else check alpha value
        if kwargs.get('doppler_hwhm_yn') is not None:
            self.doppler_hwhm_yn = kwargs['doppler_hwhm_yn']
            self.alpha_hwhm = float(_alpha) if _alpha is not None else 3.0
        elif hasattr(self, 'doppler_hwhm_yn') and 'alpha_hwhm' not in kwargs:
             # Keep file setting if alpha not overridden
             pass
        elif _alpha is None:
            self.doppler_hwhm_yn = 'N'
            self.alpha_hwhm = 3.0
        else:
            self.doppler_hwhm_yn = 'Y'
            self.alpha_hwhm = float(_alpha)
        
        # Lorentzian HWHM
        _gamma = kwargs.get('gamma_hwhm', getattr(self, 'gamma_hwhm', 0.5))
        if kwargs.get('lorentzian_hwhm_yn') is not None:
            self.lorentzian_hwhm_yn = kwargs['lorentzian_hwhm_yn']
            self.gamma_hwhm = float(_gamma) if _gamma is not None else 0.5
        elif hasattr(self, 'lorentzian_hwhm_yn') and 'gamma_hwhm' not in kwargs:
             pass
        elif _gamma is None:
            self.lorentzian_hwhm_yn = 'N'
            self.gamma_hwhm = 0.5
        else:
            self.lorentzian_hwhm_yn = 'Y'
            self.gamma_hwhm = float(_gamma)
            
        self.alpha_hwhm_colid = kwargs.get('alpha_hwhm_colid', getattr(self, 'alpha_hwhm_colid', None))
        self.gamma_hwhm_colid = kwargs.get('gamma_hwhm_colid', getattr(self, 'gamma_hwhm_colid', None))
        
        # Database-specific defaults (will be overridden by _resolve_metadata)
        if not hasattr(self, 'species_id') or self.species_id is None:
            self.species_id = 0
        self.species_main_id = getattr(self, 'species_main_id', 0)
        self.species_sub_id = getattr(self, 'species_sub_id', 0)
        self.abundance = getattr(self, 'abundance', 1.0)
        self.mass = getattr(self, 'mass', 1.0)
        self.states_col = getattr(self, 'states_col', [])
        self.states_fmt = getattr(self, 'states_fmt', [])
        self.check_uncertainty = getattr(self, 'check_uncertainty', False)
        self.check_lifetime = getattr(self, 'check_lifetime', False)
        self.check_gfactor = getattr(self, 'check_gfactor', False)
        self.check_predissoc = getattr(self, 'check_predissoc', False)
        
        # Plotting defaults
        # Oscillator Strength
        _plot_os_yn = kwargs.get('plot_oscillator_strength_yn', getattr(self, 'plot_oscillator_strength_yn', None))
        if _plot_os_yn is not None:
            self.plot_oscillator_strength_yn = _plot_os_yn
        elif 'plot_oscillator_strength' in kwargs:
            self.plot_oscillator_strength_yn = 'Y' if kwargs['plot_oscillator_strength'] else 'N'
        elif not hasattr(self, 'plot_oscillator_strength_yn'):
            self.plot_oscillator_strength_yn = 'N'

        self.plot_oscillator_strength_method = kwargs.get('plot_oscillator_strength_method', getattr(self, 'plot_oscillator_strength_method', 'log')).upper()
        self.plot_oscillator_strength_wn_wl = kwargs.get('plot_oscillator_strength_wn_wl', getattr(self, 'plot_oscillator_strength_wn_wl', 'WN')).upper()
        self.plot_oscillator_strength_unit = kwargs.get('plot_oscillator_strength_unit', getattr(self, 'plot_oscillator_strength_unit', 'cm-1')).lower()
        self.limit_yaxis_os = kwargs.get('limit_yaxis_os', getattr(self, 'limit_yaxis_os', 1e-30))
        
        # Stick Spectra
        _plot_ss_yn = kwargs.get('plot_stick_spectra_yn', getattr(self, 'plot_stick_spectra_yn', None))
        if _plot_ss_yn is not None:
            self.plot_stick_spectra_yn = _plot_ss_yn
        elif 'plot_stick_spectra' in kwargs:
            self.plot_stick_spectra_yn = 'Y' if kwargs['plot_stick_spectra'] else 'N'
        elif not hasattr(self, 'plot_stick_spectra_yn'):
            self.plot_stick_spectra_yn = 'N'
            
        self.plot_stick_spectra_method = kwargs.get('plot_stick_spectra_method', getattr(self, 'plot_stick_spectra_method', 'log')).upper()
        self.plot_stick_spectra_wn_wl = kwargs.get('plot_stick_spectra_wn_wl', getattr(self, 'plot_stick_spectra_wn_wl', 'WN')).upper()
        self.plot_stick_spectra_unit = kwargs.get('plot_stick_spectra_unit', getattr(self, 'plot_stick_spectra_unit', 'cm-1')).lower()
        self.limit_yaxis_stick_spectra = kwargs.get('limit_yaxis_stick_spectra', getattr(self, 'limit_yaxis_stick_spectra', 1e-30))
        
        # Cross Section
        _plot_xs_yn = kwargs.get('plot_cross_section_yn', getattr(self, 'plot_cross_section_yn', None))
        if _plot_xs_yn is not None:
            self.plot_cross_section_yn = _plot_xs_yn
        elif 'plot_cross_section' in kwargs:
            self.plot_cross_section_yn = 'Y' if kwargs['plot_cross_section'] else 'N'
        elif not hasattr(self, 'plot_cross_section_yn'):
            self.plot_cross_section_yn = 'N'

        self.plot_cross_section_method = kwargs.get('plot_cross_section_method', getattr(self, 'plot_cross_section_method', 'log')).upper()
        self.plot_cross_section_wn_wl = kwargs.get('plot_cross_section_wn_wl', getattr(self, 'plot_cross_section_wn_wl', 'WN')).upper()
        self.plot_cross_section_unit = kwargs.get('plot_cross_section_unit', getattr(self, 'plot_cross_section_unit', 'cm-1')).lower()
        self.limit_yaxis_xsec = kwargs.get('limit_yaxis_xsec', getattr(self, 'limit_yaxis_xsec', 1e-30))
    
    def _resolve_metadata(self):
        """Resolve database-specific metadata from definition files."""
        from pyexocross.api._resolve import resolve_database_metadata
        # Preserve already-parsed IDs (from inp_para) when available.
        current_species_main_id = getattr(self, 'species_main_id', 0)
        current_species_sub_id = getattr(self, 'species_sub_id', 0)
        meta = resolve_database_metadata(
            database=self.database,
            read_path=self.read_path,
            data_info=self.data_info,
            species_id=self.species_id,
        )
        self.states_col = meta['states_col']
        self.states_fmt = meta['states_fmt']
        self.check_uncertainty = meta['check_uncertainty']
        self.check_lifetime = meta['check_lifetime']
        self.check_gfactor = meta['check_gfactor']
        self.check_predissoc = meta['check_predissoc']
        resolved_species_main_id = meta['species_main_id']
        resolved_species_sub_id = meta['species_sub_id']
        # For ExoMol/ExoAtom, keep non-zero IDs that were already parsed from
        # inp_para(). If resolver returns zeros, do not overwrite valid values.
        if self.database in ('ExoMol', 'ExoAtom'):
            if current_species_main_id and current_species_sub_id and not (resolved_species_main_id or resolved_species_sub_id):
                self.species_main_id = current_species_main_id
                self.species_sub_id = current_species_sub_id
            else:
                self.species_main_id = resolved_species_main_id
                self.species_sub_id = resolved_species_sub_id
        else:
            self.species_main_id = resolved_species_main_id
            self.species_sub_id = resolved_species_sub_id
        self.abundance = meta['abundance']
        self.mass = meta['mass']
    
    def _set_attributes(self, params):
        """Set attributes from parsed parameters tuple."""
        # Unpack all parameters (matching inp_para return structure)
        (self.database, self.data_info, self.read_path, self.save_path, self.logs_path,
         self.conversion, self.partition_functions, self.specific_heats, self.cooling_functions,
         self.lifetimes, self.oscillator_strengths, self.stick_spectra, self.cross_sections,
         self.ncputrans, self.ncpufiles, self.chunk_size,
         self.conversion_format, self.conversion_min_freq, self.conversion_max_freq,
         self.conversion_unc, self.conversion_threshold,
         self.global_qn_label_list, self.global_qn_format_list,
         self.local_qn_label_list, self.local_qn_format_list,
         self.ntemp, self.tmax, self.compress_yn, self.gf_or_f,
         self.broadeners, self.ratios, self.T_list, self.P_list,
         self.wn_wl, self.wn_wl_unit, self.min_wnl, self.max_wnl,
         self.min_wn, self.max_wn, self.N_point, self.bin_size, self.wn_grid,
         self.predissoc_yn, self.photo, self.cutoff, self.threshold, self.unc_filter,
         self.nlte_method, self.nlte_path, self.lte_nlte,
         self.qnslabel_list, self.qnsformat_list, self.qns_label, self.qns_value,
         self.qns_format, self.qns_filter,
         self.doppler_hwhm_yn, self.lorentzian_hwhm_yn,
         self.alpha_hwhm, self.gamma_hwhm, self.alpha_hwhm_colid, self.gamma_hwhm_colid,
         self.abs_emi, self.profile,
         self.species_id, self.species_main_id, self.species_sub_id, self.abundance, self.mass,
         self.states_col, self.states_fmt,
         self.check_uncertainty, self.check_lifetime, self.check_gfactor, self.check_predissoc,
         self.plot_oscillator_strength_yn, self.plot_oscillator_strength_method,
         self.plot_oscillator_strength_wn_wl, self.plot_oscillator_strength_unit, self.limit_yaxis_os,
         self.plot_stick_spectra_yn, self.plot_stick_spectra_method,
         self.plot_stick_spectra_wn_wl, self.plot_stick_spectra_unit, self.limit_yaxis_stick_spectra,
         self.tvib_list, self.trot_list, self.vib_label, self.rot_label,
         self.plot_cross_section_yn, self.plot_cross_section_method,
         self.plot_cross_section_wn_wl, self.plot_cross_section_unit, self.limit_yaxis_xsec) = params
    
    def to_globals(self):
        """
        Set global variables from configuration.
        
        This is used to make configuration available to legacy code that uses
        global variables. Sets variables in the core module's globals.
        """
        import sys
        import pyexocross.core as core_module
        
        # Map all attributes to legacy global variable names
        globals_map = {
            'database': 'database',
            'data_info': 'data_info',
            'read_path': 'read_path',
            'save_path': 'save_path',
            'logs_path': 'logs_path',
            'conversion': 'Conversion',
            'partition_functions': 'PartitionFunctions',
            'specific_heats': 'SpecificHeats',
            'lifetimes': 'Lifetimes',
            'cooling_functions': 'CoolingFunctions',
            'oscillator_strengths': 'OscillatorStrengths',
            'stick_spectra': 'StickSpectra',
            'cross_sections': 'CrossSections',
            'ncputrans': 'ncputrans',
            'ncpufiles': 'ncpufiles',
            'chunk_size': 'chunk_size',
            'conversion_format': 'ConversionFormat',
            'conversion_min_freq': 'ConversionMinFreq',
            'conversion_max_freq': 'ConversionMaxFreq',
            'conversion_unc': 'ConversionUnc',
            'conversion_threshold': 'ConversionThreshold',
            'global_qn_label_list': 'GlobalQNLabel_list',
            'global_qn_format_list': 'GlobalQNFormat_list',
            'local_qn_label_list': 'LocalQNLabel_list',
            'local_qn_format_list': 'LocalQNFormat_list',
            'ntemp': 'Ntemp',
            'tmax': 'Tmax',
            'compress_yn': 'CompressYN',
            'gf_or_f': 'gfORf',
            'broadeners': 'broadeners',
            'ratios': 'ratios',
            'T_list': 'T_list',
            'P_list': 'P_list',
            'wn_wl': 'wn_wl',
            'wn_wl_unit': 'wn_wl_unit',
            'min_wnl': 'min_wnl',
            'max_wnl': 'max_wnl',
            'min_wn': 'min_wn',
            'max_wn': 'max_wn',
            'N_point': 'N_point',
            'bin_size': 'bin_size',
            'wn_grid': 'wn_grid',
            'predissoc_yn': 'predissocYN',
            'photo': 'photo',
            'cutoff': 'cutoff',
            'threshold': 'threshold',
            'unc_filter': 'UncFilter',
            'nlte_method': 'NLTEMethod',
            'nlte_path': 'NLTEPath',
            'lte_nlte': 'LTE_NLTE',
            'qnslabel_list': 'QNslabel_list',
            'qnsformat_list': 'QNsformat_list',
            'qns_label': 'QNs_label',
            'qns_value': 'QNs_value',
            'qns_format': 'QNs_format',
            'qns_filter': 'QNsFilter',
            'doppler_hwhm_yn': 'DopplerHWHMYN',
            'lorentzian_hwhm_yn': 'LorentzianHWHMYN',
            'alpha_hwhm': 'alpha_HWHM',
            'gamma_hwhm': 'gamma_HWHM',
            'alpha_hwhm_colid': 'alpha_hwhm_colid',
            'gamma_hwhm_colid': 'gamma_hwhm_colid',
            'abs_emi': 'abs_emi',
            'profile': 'profile',
            'species_id': 'species_id',
            'species_main_id': 'species_main_id',
            'species_sub_id': 'species_sub_id',
            'abundance': 'abundance',
            'mass': 'mass',
            'states_col': 'states_col',
            'states_fmt': 'states_fmt',
            'check_uncertainty': 'check_uncertainty',
            'check_lifetime': 'check_lifetime',
            'check_gfactor': 'check_gfactor',
            'check_predissoc': 'check_predissoc',
            'plot_oscillator_strength_yn': 'PlotOscillatorStrengthYN',
            'plot_oscillator_strength_method': 'PlotOscillatorStrengthMethod',
            'plot_oscillator_strength_wn_wl': 'PlotOscillatorStrengthWnWl',
            'plot_oscillator_strength_unit': 'PlotOscillatorStrengthUnit',
            'limit_yaxis_os': 'limitYaxisOS',
            'plot_stick_spectra_yn': 'PlotStickSpectraYN',
            'plot_stick_spectra_method': 'PlotStickSpectraMethod',
            'plot_stick_spectra_wn_wl': 'PlotStickSpectraWnWl',
            'plot_stick_spectra_unit': 'PlotStickSpectraUnit',
            'limit_yaxis_stick_spectra': 'limitYaxisStickSpectra',
            'tvib_list': 'Tvib_list',
            'trot_list': 'Trot_list',
            'vib_label': 'vib_label',
            'rot_label': 'rot_label',
            'plot_cross_section_yn': 'PlotCrossSectionYN',
            'plot_cross_section_method': 'PlotCrossSectionMethod',
            'plot_cross_section_wn_wl': 'PlotCrossSectionWnWl',
            'plot_cross_section_unit': 'PlotCrossSectionUnit',
            'limit_yaxis_xsec': 'limitYaxisXsec',
        }
        
        # Set variables in core module's globals
        for attr_name, legacy_name in globals_map.items():
            if hasattr(self, attr_name):
                value = getattr(self, attr_name)
                setattr(core_module, legacy_name, value)
                # Also set in sys.modules['__main__'] if it exists
                if '__main__' in sys.modules:
                    setattr(sys.modules['__main__'], legacy_name, value)
        
        # Also update src.base.input module-level attributes that other
        # modules access dynamically (e.g. hitran_qn.py).
        try:
            import src.base.input as _input_mod
            _qn_attrs = {
                'global_qn_label_list': 'GlobalQNLabel_list',
                'global_qn_format_list': 'GlobalQNFormat_list',
                'local_qn_label_list': 'LocalQNLabel_list',
                'local_qn_format_list': 'LocalQNFormat_list',
            }
            for attr_name, mod_name in _qn_attrs.items():
                if hasattr(self, attr_name):
                    setattr(_input_mod, mod_name, getattr(self, attr_name))
        except ImportError:
            pass
        
        # Refresh module-level variables in source modules that cache
        # get_config() results at import time.  Without this, a second
        # pyx.run() call would use stale values for database, profile,
        # DopplerHWHMYN, LorentzianHWHMYN, predissocYN, and check_predissoc.
        _stale_var_map = {
            'DopplerHWHMYN': self.doppler_hwhm_yn,
            'LorentzianHWHMYN': self.lorentzian_hwhm_yn,
            'database': self.database,
            'predissocYN': self.predissoc_yn,
            'check_predissoc': self.check_predissoc,
            'profile': self.profile,
        }
        _stale_modules = [
            'src.calculation.calcualte_line_profile',
            'calculation.calcualte_line_profile',
            'src.save.hitran.hitran_cross_section',
            'src.save.hitran.hitran_stick_spectra_cross_section',
        ]
        for _mod_name in _stale_modules:
            if _mod_name in sys.modules:
                _mod = sys.modules[_mod_name]
                for _var, _val in _stale_var_map.items():
                    setattr(_mod, _var, _val)
    
    def _to_legacy_name(self, key):
        """Convert config attribute name to legacy global variable name."""
        name_map = {
            'conversion': 'Conversion',
            'partition_functions': 'PartitionFunctions',
            'specific_heats': 'SpecificHeats',
            'lifetimes': 'Lifetimes',
            'cooling_functions': 'CoolingFunctions',
            'oscillator_strengths': 'OscillatorStrengths',
            'stick_spectra': 'StickSpectra',
            'cross_sections': 'CrossSections',
            'ncputrans': 'ncputrans',
            'ncpufiles': 'ncpufiles',
            'chunk_size': 'chunk_size',
            'conversion_format': 'ConversionFormat',
            'conversion_min_freq': 'ConversionMinFreq',
            'conversion_max_freq': 'ConversionMaxFreq',
            'conversion_unc': 'ConversionUnc',
            'conversion_threshold': 'ConversionThreshold',
            'global_qn_label_list': 'GlobalQNLabel_list',
            'global_qn_format_list': 'GlobalQNFormat_list',
            'local_qn_label_list': 'LocalQNLabel_list',
            'local_qn_format_list': 'LocalQNFormat_list',
            'ntemp': 'Ntemp',
            'tmax': 'Tmax',
            'compress_yn': 'CompressYN',
            'gf_or_f': 'gfORf',
            'broadeners': 'broadeners',
            'ratios': 'ratios',
            'T_list': 'T_list',
            'P_list': 'P_list',
            'wn_wl': 'wn_wl',
            'wn_wl_unit': 'wn_wl_unit',
            'min_wnl': 'min_wnl',
            'max_wnl': 'max_wnl',
            'min_wn': 'min_wn',
            'max_wn': 'max_wn',
            'N_point': 'N_point',
            'bin_size': 'bin_size',
            'wn_grid': 'wn_grid',
            'predissoc_yn': 'predissocYN',
            'photo': 'photo',
            'cutoff': 'cutoff',
            'threshold': 'threshold',
            'unc_filter': 'UncFilter',
            'nlte_method': 'NLTEMethod',
            'nlte_path': 'NLTEPath',
            'lte_nlte': 'LTE_NLTE',
            'qnslabel_list': 'QNslabel_list',
            'qnsformat_list': 'QNsformat_list',
            'qns_label': 'QNs_label',
            'qns_value': 'QNs_value',
            'qns_format': 'QNs_format',
            'qns_filter': 'QNsFilter',
            'doppler_hwhm_yn': 'DopplerHWHMYN',
            'lorentzian_hwhm_yn': 'LorentzianHWHMYN',
            'alpha_hwhm': 'alpha_HWHM',
            'gamma_hwhm': 'gamma_HWHM',
            'alpha_hwhm_colid': 'alpha_hwhm_colid',
            'gamma_hwhm_colid': 'gamma_hwhm_colid',
            'abs_emi': 'abs_emi',
            'profile': 'profile',
            'species_id': 'species_id',
            'species_main_id': 'species_main_id',
            'species_sub_id': 'species_sub_id',
            'abundance': 'abundance',
            'mass': 'mass',
            'states_col': 'states_col',
            'states_fmt': 'states_fmt',
            'check_uncertainty': 'check_uncertainty',
            'check_lifetime': 'check_lifetime',
            'check_gfactor': 'check_gfactor',
            'check_predissoc': 'check_predissoc',
            'plot_oscillator_strength_yn': 'PlotOscillatorStrengthYN',
            'plot_oscillator_strength_method': 'PlotOscillatorStrengthMethod',
            'plot_oscillator_strength_wn_wl': 'PlotOscillatorStrengthWnWl',
            'plot_oscillator_strength_unit': 'PlotOscillatorStrengthUnit',
            'limit_yaxis_os': 'limitYaxisOS',
            'plot_stick_spectra_yn': 'PlotStickSpectraYN',
            'plot_stick_spectra_method': 'PlotStickSpectraMethod',
            'plot_stick_spectra_wn_wl': 'PlotStickSpectraWnWl',
            'plot_stick_spectra_unit': 'PlotStickSpectraUnit',
            'limit_yaxis_stick_spectra': 'limitYaxisStickSpectra',
            'tvib_list': 'Tvib_list',
            'trot_list': 'Trot_list',
            'vib_label': 'vib_label',
            'rot_label': 'rot_label',
            'plot_cross_section_yn': 'PlotCrossSectionYN',
            'plot_cross_section_method': 'PlotCrossSectionMethod',
            'plot_cross_section_wn_wl': 'PlotCrossSectionWnWl',
            'plot_cross_section_unit': 'PlotCrossSectionUnit',
            'limit_yaxis_xsec': 'limitYaxisXsec',
        }
        return name_map.get(key, key)

