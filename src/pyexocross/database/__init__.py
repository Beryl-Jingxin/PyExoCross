"""
Database loading modules for PyExoCross.

This package provides functions for loading data from:
- ExoMol database files
- HITRAN database files
"""
from .load_exomol import (
    read_all_states,
    read_part_states,
    get_transfiles,
    get_part_transfiles,
    read_exomol_pf,
    read_broad,
    extract_broad,
)
from .load_exomolhr import (
    resolve_exomolhr_data_info,
    resolve_exomolhr_filepaths,
    exomolhr_metadata,
    read_exomolhr_df,
    read_exomolhr_pf,
    process_exomolhr_linelist,
    process_exomolhr_linelist_Q,
)
from .load_hitran import (
    read_parfile,
    read_hitran_parfile,
    read_hitran_pf,
    process_hitran_linelist,
    process_hitran_linelist_Q,
)

__all__ = [
    # ExoMol
    'read_all_states',
    'read_part_states',
    'get_transfiles',
    'get_part_transfiles',
    'read_exomol_pf',
    'read_broad',
    'extract_broad',
    # ExoMolHR
    'resolve_exomolhr_data_info',
    'resolve_exomolhr_filepaths',
    'exomolhr_metadata',
    'read_exomolhr_df',
    'read_exomolhr_pf',
    'process_exomolhr_linelist',
    'process_exomolhr_linelist_Q',
    # HITRAN
    'read_parfile',
    'read_hitran_parfile',
    'read_hitran_pf',
    'process_hitran_linelist',
    'process_hitran_linelist_Q',
]
