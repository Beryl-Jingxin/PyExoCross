"""
Format conversion modules for PyExoCross.

This package provides functions for converting between:
- ExoMol and HITRAN formats
"""
from .exomol_to_hitran import (
    read_unc_states,
    convert_QNValues_exomol2hitran,
    linelist_ExoMol2HITRAN,
    broadener_ExoMol2HITRAN,
    convert_QNFormat_exomol2hitran,
    error_code,
    convert_exomol2hitran_linelist,
    process_exomol2hitran_linelist,
    conversion_exomol2hitran,
)
from .hitran_to_exomol import (
    convert_QNValues_hitran2exomol,
    convert_hitran2StatesTrans,
    convert_hitran2broad,
    conversion_states,
    conversion_trans,
    conversion_broad,
    conversion_hitran2exomol,
)

__all__ = [
    # ExoMol to HITRAN
    'read_unc_states',
    'convert_QNValues_exomol2hitran',
    'linelist_ExoMol2HITRAN',
    'broadener_ExoMol2HITRAN',
    'convert_QNFormat_exomol2hitran',
    'error_code',
    'convert_exomol2hitran_linelist',
    'process_exomol2hitran_linelist',
    'conversion_exomol2hitran',
    # HITRAN to ExoMol
    'convert_QNValues_hitran2exomol',
    'convert_hitran2StatesTrans',
    'convert_hitran2broad',
    'conversion_states',
    'conversion_trans',
    'conversion_broad',
    'conversion_hitran2exomol',
]

