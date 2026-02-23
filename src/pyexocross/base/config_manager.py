"""
Configuration Manager for PyExoCross.

Provides a singleton pattern to cache parsed configuration parameters,
avoiding repeated parsing of the same input file.
"""
import os
from typing import Optional, Dict, Tuple


class ConfigManager:
    """
    Singleton configuration manager that caches parsed input parameters.
    
    This class ensures that each input file is only parsed once, and
    subsequent requests for the same file return the cached configuration.
    
    Examples
    --------
    >>> # First call - parses the file
    >>> params1 = ConfigManager.get_config('.input/MgH.inp')
    >>> 
    >>> # Second call - returns cached result (no parsing)
    >>> params2 = ConfigManager.get_config('.input/MgH.inp')
    >>> assert params1 is params2  # Same object
    """
    
    _instance = None
    _cache: Dict[str, Tuple] = {}
    # Keep track of the most recently loaded configuration parameters so that
    # utility modules (line profiles, save modules, etc.) can access a
    # "current" config without needing the filepath again.
    _last_params: Optional[Tuple] = None
    # Also store a Config object directly (for kwargs-based initialization).
    _last_config = None
    
    def __new__(cls):
        """Ensure only one instance exists (singleton pattern)."""
        if cls._instance is None:
            cls._instance = super(ConfigManager, cls).__new__(cls)
        return cls._instance
    
    @classmethod
    def get_config(cls, inp_filepath: str, force_reload: bool = False) -> Tuple:
        """
        Get configuration parameters for an input file.
        
        Parameters
        ----------
        inp_filepath : str
            Path to the input configuration file (.inp)
        force_reload : bool, optional
            If True, force re-parsing even if cached. Default is False.
        
        Returns
        -------
        tuple
            Parsed configuration parameters (same format as inp_para)
        
        Examples
        --------
        >>> from pyexocross.base.config_manager import ConfigManager
        >>> params = ConfigManager.get_config('.input/MgH.inp')
        >>> database, data_info, read_path, save_path = params[0:4]
        """
        # Normalize file path for caching
        abs_path = os.path.abspath(inp_filepath)
        
        # Check cache first
        if not force_reload and abs_path in cls._cache:
            params = cls._cache[abs_path]
            cls._last_params = params
            return params
        
        # Import here to avoid circular imports
        from .input import inp_para
        
        # Parse and cache
        params = inp_para(inp_filepath)
        cls._cache[abs_path] = params
        cls._last_params = params
        return params
    
    @classmethod
    def clear_cache(cls, inp_filepath: Optional[str] = None):
        """
        Clear configuration cache.
        
        Parameters
        ----------
        inp_filepath : str, optional
            If provided, clear only this file's cache.
            If None, clear all cached configurations.
        
        Examples
        --------
        >>> ConfigManager.clear_cache('.input/MgH.inp')  # Clear specific file
        >>> ConfigManager.clear_cache()  # Clear all cache
        """
        if inp_filepath is None:
            cls._cache.clear()
        else:
            abs_path = os.path.abspath(inp_filepath)
            cls._cache.pop(abs_path, None)
    
    @classmethod
    def is_cached(cls, inp_filepath: str) -> bool:
        """
        Check if a configuration is cached.
        
        Parameters
        ----------
        inp_filepath : str
            Path to the input configuration file
        
        Returns
        -------
        bool
            True if configuration is cached, False otherwise
        """
        abs_path = os.path.abspath(inp_filepath)
        return abs_path in cls._cache
    
    @classmethod
    def get_cache_size(cls) -> int:
        """
        Get the number of cached configurations.
        
        Returns
        -------
        int
            Number of cached configurations
        """
        return len(cls._cache)


# Convenience function for easy access
def get_config(inp_filepath: Optional[str] = None, force_reload: bool = False) -> Tuple:
    """
    Convenience function to get configuration parameters.
    
    This is a wrapper around ConfigManager.get_config() for easier access.
    
    Parameters
    ----------
    inp_filepath : str, optional
        Path to the input configuration file (.inp). If omitted or None,
        the most recently loaded configuration parameters are returned.
    force_reload : bool, optional
        If True, force re-parsing even if cached. Default is False.
    
    Returns
    -------
    tuple
        Parsed configuration parameters
    
    Examples
    --------
    >>> from pyexocross.base.config_manager import get_config
    >>> params = get_config('.input/MgH.inp')
    >>> database, data_info, read_path, save_path = params[0:4]
    """
    # If no filepath is provided, return the most recently loaded params.
    if inp_filepath is None:
        # Try Config object first (kwargs-based initialization)
        if ConfigManager._last_config is not None:
            cfg = ConfigManager._last_config
            return (
                cfg.doppler_hwhm_yn,
                cfg.lorentzian_hwhm_yn,
                cfg.database,
                cfg.predissoc_yn,
                cfg.check_predissoc,
                cfg.profile,
            )
        if ConfigManager._last_params is not None:
            # For backward-compatibility with legacy modules that expect
            # a *short* tuple
            # (DopplerHWHMYN, LorentzianHWHMYN, database, predissocYN,
            #  check_predissoc, profile),
            # reconstruct these values from the last full parameter set by
            # using Config._set_attributes.
            try:
                from pyexocross.config import Config as _Config
            except Exception:  # pragma: no cover - very unlikely import issue
                # Fallback: return full params if Config cannot be imported
                return ConfigManager._last_params
            cfg = _Config.__new__(_Config)
            cfg._set_attributes(ConfigManager._last_params)
            return (
                cfg.doppler_hwhm_yn,
                cfg.lorentzian_hwhm_yn,
                cfg.database,
                cfg.predissoc_yn,
                cfg.check_predissoc,
                cfg.profile,
            )
        raise ValueError(
            "No configuration has been loaded yet. "
            "Please create a Config(inp_filepath=...) first."
        )
    return ConfigManager.get_config(inp_filepath, force_reload)

