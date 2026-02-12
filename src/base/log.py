from __future__ import annotations

import os
import sys
import datetime
import atexit
import pandas as pd
from tqdm import tqdm
from .utils import Timer, ensure_dir

_ORIGINAL_STDOUT = sys.stdout
_ORIGINAL_STDERR = sys.stderr
_LOG_FILE_HANDLE = None

# Ensure print writes are flushed immediately for real-time nohup logs
try:
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)
except Exception:
    pass

class TeeStream:
    """
    Stream wrapper that writes to multiple output streams simultaneously.
    
    Useful for duplicating output to both console and log file.
    """
    
    def __init__(self, *streams):
        """
        Initialize TeeStream with multiple output streams.

        Parameters
        ----------
        *streams : file-like objects
            One or more file-like objects to write to simultaneously.
        """
        self.streams = streams

    def write(self, data):
        """
        Write data to all registered streams.

        Parameters
        ----------
        data : str or bytes
            Data to write to all streams.
        """
        for stream in self.streams:
            stream.write(data)
        self.flush()

    def flush(self):
        """
        Flush all registered streams.
        """
        for stream in self.streams:
            stream.flush()

    def isatty(self):
        """
        Check if any of the streams is a TTY (terminal).

        Returns
        -------
        bool
            True if any stream is a TTY, False otherwise.
        """
        return any(getattr(stream, 'isatty', lambda: False)() for stream in self.streams)

def _close_log_file():
    """
    Close the log file and restore original stdout/stderr streams.

    This function is registered with atexit to ensure proper cleanup.
    Restores original streams before closing to prevent write errors.
    """
    global _LOG_FILE_HANDLE
    try:
        # Restore stdout and stderr before closing the log file
        # This prevents errors when TeeStream tries to write to a closed file handle
        if sys.stdout is not _ORIGINAL_STDOUT:
            sys.stdout = _ORIGINAL_STDOUT
        if sys.stderr is not _ORIGINAL_STDERR:
            sys.stderr = _ORIGINAL_STDERR
        
        # Now safely close the log file
        if _LOG_FILE_HANDLE is not None:
            try:
                _LOG_FILE_HANDLE.flush()
            except (ValueError, OSError):
                pass  # File may already be closed or flushed
            try:
                _LOG_FILE_HANDLE.close()
            except (ValueError, OSError):
                pass  # File may already be closed
            _LOG_FILE_HANDLE = None
    except Exception:
        # Silently ignore any exceptions during cleanup
        # This prevents "Exception ignored in sys.unraisablehook" messages
        pass

def setup_logging(log_file_path):
    """
    Set up logging to both console and file with date suffix.

    Creates a log file with date suffix and redirects stdout/stderr to both
    console and file using TeeStream. Registers cleanup function for proper
    file closure on exit.

    Parameters
    ----------
    log_file_path : str
        Base path for the log file. A date suffix will be added to the filename.
    """
    global _LOG_FILE_HANDLE
    log_dir = os.path.dirname(log_file_path)
    log_base = os.path.splitext(os.path.basename(log_file_path))[0]
    log_ext = os.path.splitext(log_file_path)[1] or '.log'
    date_suffix = datetime.datetime.now().strftime('%Y%m%d')
    dated_name = f'{log_base}__{date_suffix}{log_ext}'
    final_log_path = os.path.join(log_dir, dated_name)
    # Always overwrite existing log for this date instead of appending
    _LOG_FILE_HANDLE = open(final_log_path, 'w', buffering=1, encoding='utf-8')
    sys.stdout = TeeStream(_ORIGINAL_STDOUT, _LOG_FILE_HANDLE)
    sys.stderr = TeeStream(_ORIGINAL_STDERR, _LOG_FILE_HANDLE)
    atexit.register(_close_log_file)
    print(f'Logging to file: {final_log_path}')

def parse_logging_info(inp_filepath):
    """
    Parse logging path from input file before setting up logging.
    This ensures all print statements (including those in inp_para) are logged.
    """
    # Find the maximum column for all the rows.
    with open(inp_filepath, 'r') as temp_f:
        col_count = max([len([x for x in l.split(" ") if x.strip()]) for l in temp_f.readlines()])
    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1).
    column_names = [i for i in range(col_count)] 
    inp_df = pd.read_csv(inp_filepath, sep='\\s+', header = None, names=column_names, usecols=column_names)
    col0 = inp_df[0]
    
    # File path - parse logs_path
    logs_path_raw = inp_df[col0.isin(['LogFilePath'])][1].values[0]
    logs_path_raw = logs_path_raw.replace('//','/').strip()
    if logs_path_raw == '':
        raise ValueError("LogFilePath cannot be empty.")
    log_dir = os.path.dirname(logs_path_raw)
    log_name = os.path.basename(logs_path_raw)
    if log_dir == '': 
        log_dir = os.getcwd()
    ensure_dir(log_dir + '/')
    logs_path = os.path.join(log_dir, log_name)
    return logs_path

# Parse logging info first and set up logging before calling inp_para
# This ensures all print statements (including those in inp_para) are logged
# Note: This is only executed when running as a script, not when importing as a module
# To use logging, call setup_logging() explicitly with a logs_path
# if __name__ == '__main__':
#     from .input import parse_args
#     inp_filepath = parse_args()
#     logs_path = parse_logging_info(inp_filepath)
#     setup_logging(logs_path)

class _ProgressLogger:
    """
    Progress logger that displays progress bars in both interactive and non-interactive modes.
    
    Provides a tqdm-like interface that works in both terminal and log file contexts,
    automatically adapting to the output environment. Displays progress as percentage
    bars in log files and uses tqdm for interactive terminal output.
    
    Attributes
    ----------
    iterable : iterable
        The iterable object being tracked.
    total : int or None
        Total number of items. None if unknown.
    desc : str
        Description string displayed with progress bar.
    count : int
        Current number of items processed.
    next_pct : int or None
        Next percentage threshold to print (increments of 10%).
    bar_length : int
        Length of progress bar in characters.
    log_streams : list
        List of file-like objects to write progress updates to.
    """
    
    def __init__(self, iterable, total=None, desc=''):
        """
        Initialize progress logger.

        Parameters
        ----------
        iterable : iterable
            Iterable object to track progress for.
        total : int, optional
            Total number of items. If None, progress percentage cannot be calculated.
        desc : str, optional
            Description string to display with progress bar.
        """
        self.iterable = iterable
        self.total = total
        self.desc = desc
        self.count = 0
        self.next_pct = 0 if total else None
        self.bar_length = 10
        self.log_streams = []
        if _LOG_FILE_HANDLE:
            self.log_streams.append(_LOG_FILE_HANDLE)
        if not _ORIGINAL_STDOUT.isatty():
            self.log_streams.append(sys.stdout)
        # remove duplicates while preserving order
        seen = []
        for stream in self.log_streams:
            if stream not in seen:
                seen.append(stream)
        self.log_streams = seen
        if self.total and self.next_pct == 0 and self.log_streams:
            self._print_pct(0)
            self.next_pct = 10

    def _print_pct(self, pct):
        """
        Print progress percentage to log streams.

        Formats and writes a progress bar showing the current percentage
        to all registered log streams.

        Parameters
        ----------
        pct : int
            Progress percentage (0-100) to display.
        """
        prefix = f'{self.desc} ' if self.desc else ''
        filled = min(self.bar_length, pct // 10)
        bar = '█' * filled + ' ' * (self.bar_length - filled)
        line = f'{prefix}[{bar}] {pct:3d}%\n'
        for stream in self.log_streams:
            stream.write(line)
            stream.flush()

    def __iter__(self):
        """
        Iterate over the iterable and update progress display.

        Yields items from the iterable while tracking progress and
        printing progress updates at 10% intervals.

        Yields
        ------
        object
            Items from the iterable, one at a time.
        """
        for item in self.iterable:
            self.count += 1
            if self.total and self.log_streams:
                pct = int(self.count * 100 / self.total)
                while self.next_pct is not None and pct >= self.next_pct:
                    self._print_pct(min(self.next_pct, 100))
                    self.next_pct += 10
                    if self.next_pct > 100:
                        self.next_pct = None
            yield item
        if self.total and self.next_pct is not None and self.log_streams:
            self._print_pct(100)

def log_tqdm(iterable, *args, **kwargs):
    desc = kwargs.get('desc', '')
    total = kwargs.get('total')
    if total is None:
        try:
            total = len(iterable)
        except TypeError:
            total = None
    interactive_iterable = iterable
    if _ORIGINAL_STDOUT.isatty():
        tqdm_kwargs = dict(kwargs)
        tqdm_kwargs.setdefault('leave', False)
        tqdm_kwargs.setdefault('dynamic_ncols', True)
        tqdm_kwargs.setdefault('ascii', ' █')
        interactive_iterable = tqdm(iterable, *args, file=_ORIGINAL_STDOUT, **tqdm_kwargs)
    return _ProgressLogger(interactive_iterable, total=total, desc=desc)

def print_file_info(file_name, file_col, file_fmt):
    # Find the max width for each column for alignment
    widths = [max(len(str(c)), len(str(f))) for c, f in zip(file_col, file_fmt)]
    # Print header
    print()
    print(file_name, 'file column names  :', end=' ')
    for c, w in zip(file_col, widths):
        print(f'{c:<{w}}', end='  ')
    print()
    print(file_name, 'file column formats:', end=' ')
    for f, w in zip(file_fmt, widths):
        print(f'{f:<{w}}', end='  ')
    print('\n')

# Conversion
# Print conversion information
def print_conversion_info(ConversionMinFreq, ConversionMaxFreq, GlobalQNLabel_list, GlobalQNFormat_list, 
                          LocalQNLabel_list, LocalQNFormat_list, ConversionUnc, ConversionThreshold):
    """
    Print conversion parameters and filter information.

    Displays wavenumber range, quantum number labels/formats, and any
    applied filters (uncertainty and threshold) in a formatted manner.

    Parameters
    ----------
    ConversionMinFreq : float
        Minimum wavenumber for conversion (cm⁻¹).
    ConversionMaxFreq : float
        Maximum wavenumber for conversion (cm⁻¹).
    GlobalQNLabel_list : list of str
        List of global quantum number labels.
    GlobalQNFormat_list : list of str
        List of format strings for global quantum numbers.
    LocalQNLabel_list : list of str
        List of local quantum number labels.
    LocalQNFormat_list : list of str
        List of format strings for local quantum numbers.
    ConversionUnc : float or str
        Uncertainty filter value or 'None'.
    ConversionThreshold : float or str
        Threshold filter value or 'None'.
    """
    # Print conversion wavenumber range.
    print('{:25s} : {} {} {} {} {}'.format('Wavenumber range selected', ConversionMinFreq, 'cm⁻¹', '-', ConversionMaxFreq, 'cm⁻¹'))
    # Extract global quantum numbers, print the global quantum numbers information.
    if GlobalQNLabel_list !=[]:  
        # Find the max width for each column for alignment
        widths = [max(len(str(c)), len(str(f))) for c, f in zip(GlobalQNLabel_list, GlobalQNFormat_list)]
        # Print header
        print('Selected global quantum number labels :', end=' ')
        for c, w in zip(GlobalQNLabel_list, widths):
            print(f'{c:<{w}}', end='  ')
        print()
        print('Selected global quantum number formats:', end=' ')
        for f, w in zip(GlobalQNFormat_list, widths):
            print(f'{f:<{w}}', end='  ')
        print()
    else:
        pass
    # Extract local quantum numbers, print the local quantum numbers information.
    if LocalQNLabel_list !=[]:  
        # Find the max width for each column for alignment
        widths = [max(len(str(c)), len(str(f))) for c, f in zip(LocalQNLabel_list, LocalQNFormat_list)]
        # Print header
        print('Selected  local quantum number labels :', end=' ')
        for c, w in zip(LocalQNLabel_list, widths):
            print(f'{c:<{w}}', end='  ')
        print()
        print('Selected  local quantum number formats:', end=' ')
        for f, w in zip(LocalQNFormat_list, widths):
            print(f'{f:<{w}}', end='  ')
        print()
    else:
        pass   
    # If using filter, print the filter information.
    if ConversionUnc != 'None' or ConversionThreshold != 'None':
        print('\nApply filters')
    else:
        pass
    # If uncertainty filter is applied, print the uncertainty filter information.
    if ConversionUnc != 'None':
        print('{:25s} : {:<6}'.format('Uncertainty filter', ConversionUnc), 'cm⁻¹')
    else:
        pass
    # If threshold filter is applied, print the threshold filter information.
    if ConversionThreshold != 'None':
        print('{:25s} : {:<6}'.format('Threshold filter', ConversionThreshold), 'cm/molecule')
    else:
        pass
    print()

def print_stick_info(unc_unit, threshold_unit):
    """
    Print stick spectra calculation parameters and filter information.

    Displays applied filters (quantum numbers, uncertainty, threshold),
    temperature list, wavenumber/wavelength range, NLTE method description,
    intensity type (absorption/emission), and predissociation status.

    Parameters
    ----------
    unc_unit : str
        Unit string for uncertainty (e.g., 'cm⁻¹').
    threshold_unit : str
        Unit string for threshold (e.g., 'cm/molecule').
    """
    from pyexocross.core import (  
        QNsFilter,
        threshold,
        UncFilter,
        QNs_label,
        QNs_format,
        QNs_value,
        abs_emi,
        NLTEMethod,
        T_list,
        Tvib_list,
        Trot_list,
        wn_wl,
        wn_wl_unit,
        min_wn,
        max_wn,
        min_wnl,
        max_wnl,
        predissocYN,
    )
    print('\nStick spectra calculation information')
    if QNsFilter !=[] or threshold != 'None' or UncFilter != 'None':
        print('Apply filters')
    else:
        pass
    # If quantum number filter is applied, print the quantum number filter information.
    if QNsFilter !=[]:  
        if QNs_label != []:
            # Find the max width for each column for alignment
            widths = [max(len(str(c)), len(str(f)), len(str(v))) for c, f, v in zip(QNs_label, QNs_format, QNs_value)]
            # Print header
            print('Selected quantum number labels  :', end=' ')
            for c, w in zip(QNs_label, widths):
                print(f'{c:<{w}}', end='  ')
            print()
            print('Selected quantum number formats :', end=' ')
            for f, w in zip(QNs_format, widths):
                print(f'{f:<{w}}', end='  ')
            print()
            print('Selected quantum number values  :', end=' ')
            for v, w in zip(QNs_value, widths):
                if v == ['']:
                    v = 'All'
                print(f'{str(v):<{w}}', end='  ')
            print()
        else:
            pass
    # If uncertainty filter is applied, print the uncertainty filter information.
    if UncFilter != 'None':
        print('{:25s} : {:<6}'.format('Uncertainty filter', UncFilter), unc_unit)
    else:
        pass
    # If threshold filter is applied, print the threshold filter information.
    if threshold != 'None':
        print('{:25s} : {:<6}'.format('Threshold filter', threshold), threshold_unit)
    else:
        pass
    print()
    print('{:25s} : {}'.format('Intensity', abs_emi.replace('Ab', 'Absorption').replace('Em', 'Emission')))
    NLTE_case = (NLTEMethod.replace('L', 'LTE').replace('T', 'Non-LTE').replace('D','Non-LTE').replace('P','Non-LTE'))
    NLTE_desc = (NLTEMethod.replace('L', 'Boltzmann distribution')
                 .replace('T', 'Treanor distribution with Tvib and Trot')
                 .replace('D', 'Custom vibrational density nvib and Trot')
                 .replace('P', 'Custom rovibrational population'))
    print('{:25s} : {}'.format(NLTE_case, NLTE_desc))
    print('{:25s} : {:<6}'.format('Temperatures', str(sorted(list(set(T_list)))), 'K'))
    if NLTEMethod == 'T':
        print('{:25s} : {:<6}'.format('Vibrational temperatures', str(sorted(list(set(Tvib_list)))), 'K'))
        print('{:25s} : {:<6}'.format('Rotational temperatures', str(sorted(list(set(Trot_list)))), 'K'))
    elif NLTEMethod == 'D':
        print('{:25s} : {:<6}'.format('Rotational temperatures', str(sorted(list(set(Trot_list)))), 'K'))
    elif NLTEMethod == 'P':
        pass
    if wn_wl == 'WN':
        print('{:25s} : {} {} {} {} {}'.format('Wavenumber range selected', min_wn, unc_unit, '-', max_wn, unc_unit))
    elif wn_wl == 'WL' and wn_wl_unit == 'um':
        print('{:25s} : {} {} {} {} {}'.format('Wavelength range selected', min_wnl, 'μm', '-', max_wnl, 'μm'))
    elif wn_wl == 'WL' and wn_wl_unit == 'nm':
        print('{:25s} : {} {} {} {} {}'.format('Wavelength range selected', min_wnl, 'nm', '-', max_wnl, 'nm'))
    else:
        raise ValueError("Please type the correct wavenumber or wavelength choice 'wn' or 'wl' into the input file and give the unit of wavelength in the input file.")
    if predissocYN == 'Y':
        print('{:25s} : {}'.format('Predissociation', 'Yes'))
    else:
        print('{:25s} : {}'.format('Predissociation', 'No'))

def print_xsec_info(profile_label, cutoff, UncFilter, min_wnl, max_wnl, 
                    unc_unit, threshold_unit, broad, ratio):
    """
    Print cross-section calculation information and parameters.

    Displays filter settings, calculation parameters, line profile type,
    temperature/pressure settings, and non-LTE method information.

    Parameters
    ----------
    profile_label : str
        Line profile label (e.g., 'Voigt', 'Doppler', 'Lorentzian')
    cutoff : float or str
        Wing cutoff value or 'None'
    UncFilter : float or str
        Uncertainty filter value or 'None'
    min_wnl : float
        Minimum wavenumber/wavelength
    max_wnl : float
        Maximum wavenumber/wavelength
    unc_unit : str
        Unit for uncertainty (e.g., 'cm-1')
    threshold_unit : str
        Unit for threshold
    broad : list of str
        List of broadener names
    ratio : list of float
        List of mixing ratios for broadeners
    """
    from pyexocross.core import ( 
        QNsFilter,
        threshold,
        QNs_label,
        QNs_format,
        QNs_value,
        abs_emi,
        NLTEMethod,
        T_list,
        Tvib_list,
        Trot_list,
        P_list,
        wn_wl,
        wn_wl_unit,
        min_wn,
        max_wn,
        bin_size,
        N_point,
        predissocYN,
    )
    print('\nCross section calculation information')
    # If using filter, print the filter information.
    if QNsFilter !=[] or threshold != 'None' or UncFilter != 'None':
        print('Apply filters')
    else:
        pass
    # If quantum number filter is applied, print the quantum number filter information.
    if QNsFilter !=[]:  
        if QNs_label != []:
            # Find the max width for each column for alignment
            widths = [max(len(str(c)), len(str(f)), len(str(v))) for c, f, v in zip(QNs_label, QNs_format, QNs_value)]
            # Print header
            print('Selected quantum number labels  :', end=' ')
            for c, w in zip(QNs_label, widths):
                print(f'{c:<{w}}', end='  ')
            print()
            print('Selected quantum number formats :', end=' ')
            for f, w in zip(QNs_format, widths):
                print(f'{f:<{w}}', end='  ')
            print()
            print('Selected quantum number values  :', end=' ')
            for v, w in zip(QNs_value, widths):
                if v == ['']:
                    v = 'All'
                print(f'{str(v):<{w}}', end='  ')
            print()
        else:
            pass
    # If uncertainty filter is applied, print the uncertainty filter information.
    if UncFilter != 'None':
        print('{:25s} : {:<6}'.format('Uncertainty filter', UncFilter), unc_unit)
    else:
        pass
    # If threshold filter is applied, print the threshold filter information.
    if threshold != 'None':
        print('{:25s} : {:<6}'.format('Threshold filter', threshold), threshold_unit)
    else:
        pass
    # Print the parameters information.
    print()
    print('{:25s} : {}'.format('Line profile', profile_label+' profile'))
    print('{:25s} : {}'.format('Intensity', abs_emi.replace('Ab', 'Absorption').replace('Em', 'Emission')))
    NLTE_case = (NLTEMethod.replace('L', 'LTE').replace('T', 'Non-LTE').replace('D','Non-LTE').replace('P','Non-LTE'))
    NLTE_desc = (NLTEMethod.replace('L', 'Boltzmann distribution')
                 .replace('T', 'Treanor distribution with Tvib and Trot')
                 .replace('D', 'Custom vibrational density nvib and Trot')
                 .replace('P', 'Custom rovibrational population'))
    print('{:25s} : {}'.format(NLTE_case, NLTE_desc))
    print('{:25s} : {:<6}'.format('Temperatures', str(sorted(list(set(T_list)))), 'K'))
    if NLTEMethod == 'T':
        print('{:25s} : {:<6}'.format('Vibrational temperatures', str(sorted(list(set(Tvib_list)))), 'K'))
        print('{:25s} : {:<6}'.format('Rotational temperatures', str(sorted(list(set(Trot_list)))), 'K'))
    elif NLTEMethod == 'D':
        print('{:25s} : {:<6}'.format('Rotational temperatures', str(sorted(list(set(Trot_list)))), 'K'))
    elif NLTEMethod == 'P':
        pass
    print('{:25s} : {:<6}'.format('Pressures', str(P_list)), 'bar')
    if wn_wl == 'WN':
        print('{:25s} : {} {} {} {} {}'.format('Wavenumber range selected', min_wn, unc_unit, '-', max_wn, unc_unit))
        print('{:25s} : {:<6}'.format('Bin size', bin_size), unc_unit)
    elif wn_wl == 'WL' and wn_wl_unit == 'um':
        print('{:25s} : {} {} {} {} {}'.format('Wavelength range selected', min_wnl, 'μm', '-', max_wnl, 'μm'))
        print('{:25s} : {:<6}'.format('Bin size', bin_size), 'μm')
    elif wn_wl == 'WL' and wn_wl_unit == 'nm':
        print('{:25s} : {} {} {} {} {}'.format('Wavelength range selected', min_wnl, 'nm', '-', max_wnl, 'nm'))
        print('{:25s} : {:<6}'.format('Bin size', bin_size), 'nm')
    else:
        raise ValueError("Please type the correct wavenumber or wavelength choice 'wn' or 'wl' into the input file and give the unit of wavelength in the input file.")
    print('{:25s} : {:<6}'.format('Number of points', N_point))
    if cutoff != 'None':
        print('{:25s} : {:<6}'.format('Wing cutoff', cutoff), unc_unit)
    else:
        print('{:25s} : {:<6}'.format('Wing cutoff', 'None'), unc_unit)
    if predissocYN == 'Y':
        print('{:25s} : {}'.format('Predissociation', 'Yes'))
    else:
        print('{:25s} : {}'.format('Predissociation', 'No'))
    print('{:25s} : {}'.format('Broadeners', str(broad).replace('[','').replace(']','').replace("'",'')))
    print('{:25s} : {}\n'.format('Ratios', str(ratio).replace('[','').replace(']','')))

def print_T_Tvib_Trot_P_path_info(T, Tvib, Trot, P, abs_emi, NLTEMethod, stick_xsec_str, file_path):
    """
    Print temperature, vibrational temperature, and rotational temperature information.

    Parameters
    ----------
    T : float
        Temperature in unit of K
    Tvib : float
        Vibrational temperature in unit of K
    Trot : float
        Rotational temperature in unit of K
    """
    P_str = f', P={P} bar' if P is not None else ''
    file_path_str = f': {file_path}' if file_path is not None else ''
    if abs_emi == 'Ab':
        if NLTEMethod == 'L' or NLTEMethod == 'P':
            print(f'{stick_xsec_str} file saved for T={T} K{P_str}{file_path_str}')
        elif NLTEMethod == 'T':
            print(f'{stick_xsec_str} file saved for Tvib={Tvib} K, Trot={Trot} K{P_str}{file_path_str}')
        elif NLTEMethod == 'D':
            print(f'{stick_xsec_str} file saved for T={T} K, Trot={Trot} K{P_str}{file_path_str}')
        else:
            raise ValueError("Please choose one LTE or non-LTE method from: 'L', 'T', 'D' or 'P'.")
    elif abs_emi == 'Em':
        if NLTEMethod == 'L':
            print(f'{stick_xsec_str} file saved for T={T} K{P_str}{file_path_str}')
        elif NLTEMethod == 'T':
            print(f'{stick_xsec_str} file saved for Tvib={Tvib} K, Trot={Trot} K{P_str}{file_path_str}')
        elif NLTEMethod == 'D':
            print(f'{stick_xsec_str} file saved for Trot={Trot} K{P_str}{file_path_str}')
        elif NLTEMethod == 'P':
            print(f'{stick_xsec_str} file saved{P_str}{file_path_str}')
        else:
            raise ValueError("Please choose one LTE or non-LTE method from: 'L', 'T', 'D' or 'P'.")
    else:
        raise ValueError("Please choose one from: 'Absorption' or 'Emission'.")
    