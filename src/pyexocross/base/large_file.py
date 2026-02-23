"""
Utilities for handling large transition files.

This module provides functions for efficiently processing large transition files
that may not fit in memory, using chunked reading and parallel processing.
"""
import os
import bz2
import subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED

from .constants import (
    DEFAULT_CHUNK_SIZE,
    LARGE_TRANS_FILE_BYTES,
    MAX_LARGE_FILE_WORKERS,
    MAX_INFLIGHT_MULTIPLIER,
    LARGE_WRITE_CHUNK_ROWS,
)
from .utils import ensure_dir

def is_large_trans_file(trans_filepath):
    """
    Check if a transition file exceeds the large file size threshold.

    Parameters
    ----------
    trans_filepath : str
        Path to the transition file to check.

    Returns
    -------
    bool
        True if the file size is greater than or equal to LARGE_TRANS_FILE_BYTES,
        False otherwise or if an error occurs accessing the file.
    """
    try:
        return os.path.getsize(trans_filepath) >= LARGE_TRANS_FILE_BYTES
    except OSError:
        return False

# Decompress Large .trans.bz2 Files
def command_decompress(trans_filename):
    """
    Decompress a bzip2-compressed transition file.

    Decompresses the file to a dedicated decompressed directory if it doesn't
    already exist. Uses system bunzip2 command for decompression.

    Parameters
    ----------
    trans_filename : str
        Path to the compressed .trans.bz2 file to decompress.

    Returns
    -------
    tuple of (str, int)
        A tuple containing:
        - Path to the decompressed file
        - Flag indicating if decompression was performed (1) or file already existed (0)
    """
    # Directory where the decompressed .trans files will be saved
    trans_dir = read_path + '/'.join(data_info) + '/decompressed/'
    ensure_dir(trans_dir)
    trans_file = os.path.join(trans_dir, trans_filename.split('/')[-1].replace('.bz2', ''))
    if os.path.exists(trans_file):
        num = 0
    else:
        command = f'bunzip2 < {trans_filename} > {trans_file}'
        print('Decompressing file:', trans_filename)
        subprocess.run(command, shell=True)
        num = 1
    return(trans_file, num)

def read_trans_chunks(trans_filepath, usecols, names, chunk_sz=None):
    """
    Read a transition file in chunks for memory-efficient processing.

    Reads the file as an iterator of DataFrame chunks, automatically handling
    bzip2 compression if the file extension indicates it.

    Parameters
    ----------
    trans_filepath : str
        Path to the transition file to read.
    usecols : list of int
        Column indices to read from the file.
    names : list of str
        Column names to assign to the selected columns.
    chunk_sz : int, optional
        Number of rows per chunk. If None, uses the global chunk_size variable
        or DEFAULT_CHUNK_SIZE as fallback.

    Returns
    -------
    TextFileReader
        A pandas TextFileReader iterator that yields DataFrame chunks.
    """
    if chunk_sz is None:
        chunk_sz = globals().get('chunk_size', DEFAULT_CHUNK_SIZE)
    read_kwargs = dict(sep='\s+', header=None, usecols=usecols, names=names,
                       chunksize=chunk_sz, iterator=True, low_memory=False, encoding='utf-8')
    if trans_filepath.endswith('.bz2'):
        read_kwargs['compression'] = 'bz2'
    return pd.read_csv(trans_filepath, **read_kwargs)

def process_large_chunks(trans_reader, handler, combine_fn, zero_factory, desc,
                         max_workers=MAX_LARGE_FILE_WORKERS, max_inflight=None, reducer=None,
                         executor_class=ProcessPoolExecutor):
    """
    Process large transition file chunks in parallel with bounded memory usage.

    Processes chunks from a file reader using a process pool executor, with
    controlled parallelism to limit memory consumption. Supports both accumulation
    and collection modes for combining results.

    Parameters
    ----------
    trans_reader : iterable
        Iterator that yields DataFrame chunks to process.
    handler : callable
        Function to process each chunk. Should accept a DataFrame chunk and return a result.
    combine_fn : callable
        Function to combine multiple results into a single result (used when reducer is None).
    zero_factory : callable
        Function that returns an empty result of the appropriate type.
    desc : str
        Description string for progress bar display.
    max_workers : int, optional
        Maximum number of worker processes. Default is MAX_LARGE_FILE_WORKERS.
    max_inflight : int, optional
        Maximum number of concurrent futures. If None, calculated as
        max_workers * MAX_INFLIGHT_MULTIPLIER.
    reducer : callable, optional
        Function to reduce two results into one. If provided, uses accumulation
        mode instead of collection mode. Should accept (accumulator, new_result)
        and return combined result.

    Returns
    -------
    object
        Combined result from all chunks, or empty result from zero_factory if no chunks processed.
    """
    if max_inflight is None:
        max_inflight = max_workers * MAX_INFLIGHT_MULTIPLIER
    results = [] if reducer is None else None
    accumulator = None
    with executor_class(max_workers=max_workers) as executor:
        futures = set()
        for trans_df_chunk in tqdm(trans_reader, desc=desc):
            futures.add(executor.submit(handler, trans_df_chunk))
            if len(futures) >= max_inflight:
                done, futures = wait(futures, return_when=FIRST_COMPLETED)
                for done_future in done:
                    result = done_future.result()
                    if reducer is None:
                        results.append(result)
                    else:
                        accumulator = result if accumulator is None else reducer(accumulator, result)
        for future in as_completed(futures):
            result = future.result()
            if reducer is None:
                results.append(result)
            else:
                accumulator = result if accumulator is None else reducer(accumulator, result)
    if reducer is None:
        return combine_fn(results) if results else zero_factory()
    return accumulator if accumulator is not None else zero_factory()

def _prepare_array_writer(data):
    """
    Prepare data for chunked writing by creating a row getter function.

    Converts input data (DataFrame or array-like) into a format suitable for
    chunked writing, returning the total number of rows and a getter function.

    Parameters
    ----------
    data : pd.DataFrame or array-like
        Data to prepare for writing. Can be a DataFrame or any array-like structure.

    Returns
    -------
    tuple of (int, callable)
        A tuple containing:
        - Total number of rows in the data
        - Getter function that accepts (start, end) indices and returns a numpy array slice
    """
    if isinstance(data, pd.DataFrame):
        total_rows = len(data)
        def getter(start, end):
            return data.iloc[start:end].to_numpy()
    else:
        arr = np.asarray(data)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        total_rows = arr.shape[0]
        def getter(start, end):
            return arr[start:end]
    return total_rows, getter

def save_large_txt(filepath_or_buffer, data, fmt, header=None, comments='', chunk_rows=LARGE_WRITE_CHUNK_ROWS, mode='w'):
    """
    Save large arrays or DataFrames to text file in chunks to manage memory.

    Writes data to a file in chunks to avoid loading the entire dataset into memory.
    Supports both regular files and bzip2-compressed files. Can write to a file path
    or an open file-like object.

    Parameters
    ----------
    filepath_or_buffer : str or file-like object
        Path to output file or an open file handle. If string ends with '.bz2',
        file will be compressed with bzip2.
    data : pd.DataFrame or array-like
        Data to write. Can be a DataFrame or any array-like structure.
    fmt : str or sequence of str
        Format string(s) for numpy.savetxt formatting.
    header : str, optional
        Header string to write at the beginning of the file. Default is None.
    comments : str, optional
        Comment prefix for header line. Default is empty string.
    chunk_rows : int, optional
        Number of rows to write per chunk. Default is LARGE_WRITE_CHUNK_ROWS.
    mode : str, optional
        File open mode ('w' for write, 'a' for append). Default is 'w'.
    """
    total_rows, getter = _prepare_array_writer(data)
    if hasattr(filepath_or_buffer, 'write'):
        f = filepath_or_buffer
        close_after = False
    else:
        # Handle bz2 compressed files
        if isinstance(filepath_or_buffer, str) and filepath_or_buffer.endswith('.bz2'):
            f = bz2.open(filepath_or_buffer, mode + 't')  # 't' for text mode
        else:
            f = open(filepath_or_buffer, mode)
        close_after = True
    try:
        if header is not None and header != '':
            f.write((comments or '') + header + '\n')
        if total_rows == 0:
            np.savetxt(f, np.empty((0, 1)), fmt=fmt)
            return
        for start in range(0, total_rows, chunk_rows):
            end = min(start + chunk_rows, total_rows)
            chunk = getter(start, end)
            if chunk.size == 0:
                continue
            np.savetxt(f, chunk, fmt=fmt)
    finally:
        if close_after:
            f.close()
            