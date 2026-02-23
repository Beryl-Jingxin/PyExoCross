"""
Base utilities and common functionality for PyExoCross.

This module provides:
- Constants and physical parameters
- Utility functions (file I/O, directory management)
- Logging and progress tracking
- Input parsing
- Large file handling
"""
from .constants import (
    DEFAULT_CHUNK_SIZE,
    LARGE_TRANS_FILE_BYTES,
    MAX_LARGE_FILE_WORKERS,
    MAX_INFLIGHT_MULTIPLIER,
    LARGE_WRITE_CHUNK_ROWS,
    num_cpus,
    Tref,
    Pref,
    N_A,
    h,
    c,
    kB,
    R,
    c2,
    c2InvTref,
    PI,
    hc,
    ln22,
    sinPI,
    SqrtPI,
    Sqrtln2,
    OneminSqrtPIln2,
    Negln2,
    PI4c,
    Inv8Pic,
    hcInv4Pi,
    Inv2ln2,
    InvSqrt2,
    InvSqrtPi,
    InvSprtln2,
    InvSqrt2Pi,
    InvSqrt2ln2,
    TwoSqrt2ln2,
    Sqrtln2InvPi,
    get_doppler_constants,
    get_bin_size_constants,
)
from .utils import ensure_dir, Timer
from .log import (
    TeeStream,
    setup_logging,
    parse_logging_info,
    _ProgressLogger,
    log_tqdm,
    print_file_info,
    print_conversion_info,
    print_stick_info,
)
from .input import parse_args, parse_TP_values, inp_para
from .config_manager import ConfigManager, get_config
from .large_file import (
    is_large_trans_file,
    command_decompress,
    read_trans_chunks,
    process_large_chunks,
    _prepare_array_writer,
    save_large_txt,
)

__all__ = [
    # Constants
    'DEFAULT_CHUNK_SIZE',
    'LARGE_TRANS_FILE_BYTES',
    'MAX_LARGE_FILE_WORKERS',
    'MAX_INFLIGHT_MULTIPLIER',
    'LARGE_WRITE_CHUNK_ROWS',
    'num_cpus',
    'Tref',
    'Pref',
    'N_A',
    'h',
    'c',
    'kB',
    'R',
    'c2',
    'c2InvTref',
    'PI',
    'hc',
    'ln22',
    'sinPI',
    'SqrtPI',
    'Sqrtln2',
    'OneminSqrtPIln2',
    'Negln2',
    'PI4c',
    'Inv8Pic',
    'hcInv4Pi',
    'Inv2ln2',
    'InvSqrt2',
    'InvSqrtPi',
    'InvSprtln2',
    'InvSqrt2Pi',
    'InvSqrt2ln2',
    'TwoSqrt2ln2',
    'Sqrtln2InvPi',
    'get_doppler_constants',
    'get_bin_size_constants',
    # Utils
    'ensure_dir',
    'Timer',
    # Logging
    'TeeStream',
    'setup_logging',
    'parse_logging_info',
    '_ProgressLogger',
    'log_tqdm',
    'print_file_info',
    'print_conversion_info',
    'print_stick_info',
    # Input
    'parse_args',
    'parse_TP_values',
    'inp_para',
    # Config Manager
    'ConfigManager',
    'get_config',
    # Large file
    'is_large_trans_file',
    'command_decompress',
    'read_trans_chunks',
    'process_large_chunks',
    '_prepare_array_writer',
    'save_large_txt',
]

