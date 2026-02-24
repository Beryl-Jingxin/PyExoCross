"""
PyExoCross: A Python package for calculating molecular cross sections and spectra.

PyExoCross provides tools for calculating absorption/emission cross sections,
stick spectra, partition functions, and other molecular properties from ExoMol,
ExoAtom, HITRAN, and HITEMP databases.

Recommended import convention::

    import pyexocross as px

Quick start::

    # Run all functions from an .inp file (equivalent to CLI)
    px.run('/path/to/config.inp')

    # Or call individual functions with keyword arguments
    px.cross_sections(database='ExoMol', molecule='MgH', ...)

Functions
---------
- ``px.run(inp_filepath)``         -- Run all functions from .inp file
- ``px.conversion(...)``           -- Format conversion
- ``px.partition_functions(...)``  -- Partition functions
- ``px.specific_heats(...)``       -- Specific heats
- ``px.cooling_functions(...)``    -- Cooling functions
- ``px.lifetimes(...)``            -- Radiative lifetimes
- ``px.oscillator_strengths(...)`` -- Oscillator strengths
- ``px.stick_spectra(...)``        -- Stick spectra
- ``px.cross_sections(...)``       -- Cross sections
"""

__version__ = "1.0.0.dev2"

from pyexocross.api import (
    # Run from .inp file
    run,
    # Primary API (snake_case plurals)
    conversion,
    partition_functions,
    specific_heats,
    cooling_functions,
    lifetimes,
    oscillator_strengths,
    stick_spectra,
    cross_sections,
    # Legacy aliases kept for backward compatibility
    convert_exomol_to_hitran,
    convert_hitran_to_exomol,
    partition_function,
    specific_heat,
    lifetime,
    cooling_function,
    oscillator_strength,
    stick_spectra,
    cross_section,
)

__all__ = [
    # Run from .inp
    'run',
    # Primary API
    'conversion',
    'partition_functions',
    'specific_heats',
    'cooling_functions',
    'lifetimes',
    'oscillator_strengths',
    'stick_spectra',
    'cross_sections',
    # Legacy aliases
    'convert_exomol_to_hitran',
    'convert_hitran_to_exomol',
    'partition_function',
    'specific_heat',
    'lifetime',
    'cooling_function',
    'oscillator_strength',
    'stick_spectra',
    'cross_section',
]
