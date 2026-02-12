"""
Input parsing utilities for PyExoCross.

This module parses the ``.inp`` configuration file and converts it into Python
parameters used by the rest of the codebase.
"""

import argparse
import datetime
import glob
import os
import re
import warnings

import numpy as np
import pandas as pd

from .constants import num_cpus
from .utils import ensure_dir


class InputWarning(UserWarning):
    pass

# Module-level defaults for quantum number labels/formats.
# These are populated at runtime by Config.to_globals() so that modules
# which import them (e.g. hitran_qn) always see the current values.
GlobalQNLabel_list = []
GlobalQNFormat_list = []
LocalQNLabel_list = []
LocalQNFormat_list = []

# The input file path
def parse_args():
    """
    Parse command-line arguments to get the input file path.

    Returns
    -------
    str
        Path to the input configuration file.
    """
    parse = argparse.ArgumentParser(description='PyExoCross Program')
    parse.add_argument('-p', '--path', type=str, metavar='', required=True, help='Input file path')
    args = parse.parse_args()
    inp_filepath = args.path
    return inp_filepath

# Process temperature and pressure values.
def parse_TP_values(value_str):
    """
    Parse temperature or pressure values from string input.

    Supports three formats:
    - Single value: "100"
    - Comma-separated list: "100,200,300"
    - Range: "100:1000:50" (start:end:step)

    Parameters
    ----------
    value_str : str
        String containing temperature/pressure values in one of the supported formats.

    Returns
    -------
    list of float
        List of parsed numeric values.

    Raises
    ------
    ValueError
        If range step is zero or range format is invalid.
    """
    # Parse temperature and pressure values, support single value, list or range.
    if ',' in value_str:
        # Multiple values separated by commas.
        return [float(x.strip()) for x in value_str.split(',')]
    elif ':' in value_str:
        # Range format: start:end:step
        parts = value_str.split(':')
        if len(parts) == 3:
            start, end, step = map(float, parts)
            if step == 0:
                raise ValueError('Range step cannot be zero.')
            return np.arange(start, end + step, step).tolist()
        else:
            raise ValueError("Range format should be 'start:end:step'")
    else:
        # Single value.
        return [float(value_str)]
    
def inp_para(inp_filepath):
    """
    Parse input configuration file and extract all program parameters.

    Reads the input file, extracts database information, file paths, function
    flags, computational parameters, filter settings, and other configuration
    options. Returns a comprehensive tuple of all parsed parameters.

    Parameters
    ----------
    inp_filepath : str
        Path to the input configuration file.

    Returns
    -------
    tuple
        A large tuple containing all parsed parameters in the following order:
        - Database and data info (database, data_info, read_path, save_path, logs_path)
        - Function flags (Conversion, PartitionFunctions, SpecificHeats, etc.)
        - Computational parameters (ncputrans, ncpufiles, chunk_size, etc.)
        - Conversion parameters (if applicable)
        - Temperature/pressure lists (T_list, P_list)
        - Wavenumber/wavelength settings
        - Filter parameters (threshold, UncFilter, QNsFilter, etc.)
        - NLTE parameters (if applicable)
        - Database-specific parameters (species_id, abundance, mass, etc.)
        - Plotting parameters (if applicable)
        See function return statement for complete list.

    Raises
    ------
    ValueError
        If required parameters are missing or invalid.
    FileNotFoundError
        If required database definition files cannot be found.
    """
    print('Date:', datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    # Find the maximum column for all the rows.
    with open(inp_filepath, 'r') as temp_f:
        col_count = max([len([x for x in l.split(" ") if x.strip()]) for l in temp_f.readlines()])
    # Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1).
    column_names = [i for i in range(col_count)] 
    inp_df = pd.read_csv(inp_filepath, sep='\\s+', header = None, names=column_names, usecols=column_names)
    col0 = inp_df[0]
    
    # Database
    database = (inp_df[col0.isin(['Database'])][1].values[0].upper()
                .replace('EXOMOL','ExoMol').replace('EXOATOM','ExoAtom'))
    
    # Basic information
    if database != 'ExoAtom':
        molecule = inp_df[col0.isin(['Molecule'])][1].values[0]
        isotopologue = inp_df[col0.isin(['Isotopologue'])][1].values[0]
        dataset = inp_df[col0.isin(['Dataset'])][1].values[0]
        data_info = [molecule, isotopologue, dataset]
    elif database == 'ExoAtom':
        atom = inp_df[col0.isin(['Atom'])][1].values[0]
        dataset = inp_df[col0.isin(['Dataset'])][1].values[0]
        data_info = [atom, dataset]
    else:
        raise ValueError("Please type the correct database choice 'ExoMol', 'ExoAtom', 'HITRAN', or 'HITEMP' into the input file.")
    try:
        species_id = int(inp_df[col0.isin(['SpeciesID'])][1].iloc[0])
    except:
        species_id = 0
    species_main_id = int(species_id/10)
    species_sub_id = int(species_id - species_main_id * 10)

    # File path
    read_path = (inp_df[col0.isin(['ReadPath'])][1].values[0] + '/').replace('//','/')
    save_path = (inp_df[col0.isin(['SavePath'])][1].values[0] + '/').replace('//','/')
    logs_path_raw = inp_df[col0.isin(['LogFilePath'])][1].values[0]
    logs_path_raw = logs_path_raw.replace('//','/').strip()
    if logs_path_raw == '':
        raise ValueError("LogFilePath cannot be empty.")
    log_dir = os.path.dirname(logs_path_raw)
    log_name = os.path.basename(logs_path_raw)
    if log_dir == '':
        log_dir = os.getcwd()
    ensure_dir(save_path)
    ensure_dir(log_dir + '/')
    logs_path = os.path.join(log_dir, log_name)
        
    # Functions 
    Conversion = int(inp_df[col0.isin(['Conversion'])][1].iloc[0])
    PartitionFunctions = int(inp_df[col0.isin(['PartitionFunctions'])][1].iloc[0])
    CoolingFunctions = int(inp_df[col0.isin(['CoolingFunctions'])][1].iloc[0])
    Lifetimes = int(inp_df[col0.isin(['Lifetimes'])][1].iloc[0])
    OscillatorStrengths = int(inp_df[col0.isin(['OscillatorStrengths'])][1].iloc[0])
    SpecificHeats = int(inp_df[col0.isin(['SpecificHeats'])][1].iloc[0])
    StickSpectra = int(inp_df[col0.isin(['StickSpectra'])][1].iloc[0])
    CrossSections = int(inp_df[col0.isin(['CrossSections'])][1].iloc[0])
    
    # Cores and chunks
    ncputrans = int(inp_df[col0.isin(['NCPUtrans'])][1].iloc[0])
    ncpufiles = int(inp_df[col0.isin(['NCPUfiles'])][1].iloc[0])
    chunk_size = int(inp_df[col0.isin(['ChunkSize'])][1].iloc[0])
    used_cpus = min(num_cpus, ncpufiles * ncputrans)
    print('Number of CPU on system:', num_cpus)
    print('Number of CPU in used  :', used_cpus, '\n')
    
    # Quantum numbers
    NeedQNs = Conversion + StickSpectra + CrossSections
    if NeedQNs != 0:
        QNslabel_list = list(inp_df[col0.isin(['QNslabel'])].iloc[0])[1:]
        QNsformat_list = list(inp_df[col0.isin(['QNsformat'])].iloc[0])[1:]
        QNslabel_list = [x for x in QNslabel_list if x == x]
        QNsformat_list = [x for x in QNsformat_list if x == x]
    else:
        QNslabel_list = []
        QNsformat_list = []  
    
    # Convert from one format to another
    if Conversion != 0:
        ConversionFormat = int(inp_df[col0.isin(['ConversionFormat'])][1].iloc[0])
        ConversionMinFreq = float(inp_df[col0.isin(['ConversionFrequncyRange'])][1].iloc[0])
        ConversionMaxFreq = float(inp_df[col0.isin(['ConversionFrequncyRange'])][2].iloc[0])
        GlobalQNLabel_list = list(inp_df[col0.isin(['GlobalQNLabel'])].iloc[0].dropna())[1:]
        GlobalQNFormat_list = list(inp_df[col0.isin(['GlobalQNFormat'])].iloc[0].dropna())[1:]
        LocalQNLabel_list = list(inp_df[col0.isin(['LocalQNLabel'])].iloc[0].dropna())[1:]
        LocalQNFormat_list = list(inp_df[col0.isin(['LocalQNFormat'])].iloc[0].dropna())[1:]
        # Uncertainty filter
        ConversionUncYN = inp_df[col0.isin(['ConvUncFilter(Y/N)'])][1].values[0].upper()[0]
        if ConversionUncYN == 'Y':
            ConversionUnc = float(inp_df[col0.isin(['ConvUncFilter(Y/N)'])][2].iloc[0])
        elif ConversionUncYN == 'N':
            ConversionUnc = 'None'
        else:
            raise ValueError("Please type the correct uncertainty filter choice 'Y' or 'N' into the input file.")  
        # Threshold filter
        ConversionThresholdYN = inp_df[col0.isin(['ConvThreshold(Y/N)'])][1].values[0].upper()[0]
        if ConversionThresholdYN == 'Y':
            ConversionThreshold = float(inp_df[col0.isin(['ConvThreshold(Y/N)'])][2].iloc[0])
        elif ConversionThresholdYN == 'N':
            ConversionThreshold = 'None'
        else:
            raise ValueError("Please type the correct threshold choice 'Y' or 'N' into the input file.")       
    else:
        ConversionFormat = 0
        ConversionMinFreq = 0
        ConversionMaxFreq = 1e20
        GlobalQNLabel_list = []
        GlobalQNFormat_list = []
        LocalQNLabel_list = []
        LocalQNFormat_list = []
        ConversionUnc = 'None'
        ConversionThreshold = 'None'
        
    # Calculate partition, cooling functions or specific heats 
    if PartitionFunctions + CoolingFunctions + SpecificHeats != 0:
        Ntemp = int(inp_df[col0.isin(['Ntemp'])][1].iloc[0])    # The number of temperature steps
        Tmax = int(inp_df[col0.isin(['Tmax'])][1].iloc[0])      # Maximal temperature in K (minimal T = 1 K )
    else:
        Ntemp = 0
        Tmax = 0  
     
    # Calculate lifetimes 
    if Lifetimes != 0:
        CompressYN = inp_df[col0.isin(['Compress(Y/N)'])][1].values[0].upper()[0]

    else:
        CompressYN = 'N'
        
    # Calculate oscillator strengths
    if OscillatorStrengths != 0:
        gfORf = inp_df[col0.isin(['gf/f'])][1].values[0].upper()
        PlotOscillatorStrengthYN = inp_df[col0.isin(['PlotOscillatorStrength(Y/N)'])][1].values[0].upper()[0]
        PlotOscillatorStrengthMethod = inp_df[col0.isin(['PlotOscillatorStrengthMethod'])][1].values[0].upper()
        PlotOscillatorStrengthWnWl = inp_df[col0.isin(['PlotOscillatorStrengthWnWl'])][1].values[0].upper()
        PlotOscillatorStrengthUnit = inp_df[col0.isin(['PlotOscillatorStrengthWnWl'])][2].values[0].lower()
        if PlotOscillatorStrengthYN == 'Y':
            _limitYaxisOS = inp_df[col0.isin(['Y-axisLimitOscillatorStrength'])][1].values[0]
            if _limitYaxisOS == '#':
                limitYaxisOS = 1e-30
            elif pd.isnull(_limitYaxisOS) == True:
                limitYaxisOS = 1e-30
            else:
                limitYaxisOS = float(_limitYaxisOS)
        else:
            limitYaxisOS = 0
    else:
        gfORf = 'GF'
        PlotOscillatorStrengthYN = 'None'
        PlotOscillatorStrengthMethod = 'None'
        PlotOscillatorStrengthWnWl = 'None'
        PlotOscillatorStrengthUnit = 'None'
        limitYaxisOS = 0
    
    # Calculate stick spectra or cross sections
    if StickSpectra + CrossSections != 0:
        LTEORNLTE = inp_df[col0.isin(['LTE/Non-LTE'])][1].values[0].upper()[0]
        # Parse multiple temperature values
        temp_str = inp_df[col0.isin(['Temperatures'])][1].values[0]
        T_list = [int(temp) for temp in parse_TP_values(temp_str)]
        wn_wl = inp_df[col0.isin(['WnWlUnit'])][1].values[0].upper()
        wn_wl_unit = inp_df[col0.isin(['WnWlUnit'])][2].values[0].lower()
        min_wnl = float(inp_df[col0.isin(['Range'])][1].iloc[0])
        max_wnl = float(inp_df[col0.isin(['Range'])][2].iloc[0])
        if wn_wl == 'WN':
            min_wn = min_wnl
            max_wn = max_wnl
        elif wn_wl == 'WL' and wn_wl_unit == 'um':
            min_wn = 1e4 / max_wnl
            try:
                max_wn = 1e4 / min_wnl
            except:
                max_wn = 1e4 / 1e-6
        elif wn_wl == 'WL' and wn_wl_unit == 'nm':
            min_wn = 1e7 / max_wnl
            try:
                max_wn = 1e7 / min_wnl
            except:
                max_wn = 1e7 / 1e-6
        else:
            raise ValueError("Please type the correct wavenumber or wavelength choice 'wn' or 'wl' and give the unit of wavelength in the input file.")
        abs_emi = inp_df[col0.isin(['Absorption/Emission'])][1].values[0].upper()[0].replace('A','Ab').replace('E','Em')
        # Uncertainty filter
        UncFilterYN = inp_df[col0.isin(['UncFilter(Y/N)'])][1].values[0].upper()[0]
        if UncFilterYN == 'Y':
            UncFilter = float(inp_df[col0.isin(['UncFilter(Y/N)'])][2].iloc[0])
        elif UncFilterYN == 'N':
            UncFilter = 'None'
        else:
            raise ValueError("Please type the correct uncertainty filter choice 'Y' or 'N' into the input file.")  
        # Threshold filter
        thresholdYN = inp_df[col0.isin(['Threshold(Y/N)'])][1].values[0].upper()[0]
        if thresholdYN == 'Y':
            threshold = float(inp_df[col0.isin(['Threshold(Y/N)'])][2].iloc[0])
        elif thresholdYN == 'N':
            threshold = 'None'
        else:
            raise ValueError("Please type the correct threshold choice 'Y' or 'N' into the input file.") 
        # Quantum number filter
        QNsFilterYN = inp_df[col0.isin(['QNsFilter(Y/N)'])][1].values[0].upper()[0]
        if QNsFilterYN == 'Y':
            QNsFilter = list(inp_df[col0.isin(['QNsFilter(Y/N)'])].iloc[0].dropna())[2:]
            QNs_label = []
            QNs_value = []
            for i in range(len(QNsFilter)):
                QNs_label.append(QNsFilter[i].split('[')[0])
                if '[' in QNsFilter[i]:
                    QNs_value.append(QNsFilter[i].split('[')[1].split(']')[0].split(';'))
                else:
                    # No bracket means select all values for this QN
                    QNs_value.append([''])
            QNs_format = [QNsformat_list[j] for j in [QNslabel_list.index(i) for i in QNs_label]]
        elif QNsFilterYN == 'N':
            QNsFilter = []
            QNs_label = []
            QNs_value = []
            QNs_format = []
        else:
            raise ValueError("Please type the correct quantum number filter choice 'Y' or 'N' into the input file.")
    else:
        LTEORNLTE = 'None'
        T_list = []
        wn_wl = 'None'
        wn_wl_unit = 'None'
        min_wnl = 0
        max_wnl = 1e20
        min_wn = 0
        max_wn = 1e20
        abs_emi = 'None'
        UncFilter = 'None'
        threshold = 'None'
        QNsFilter = []
        QNs_label = []
        QNs_value = []  
        QNs_format = []
        
    # Non-LTE
    if LTEORNLTE == 'N':
        # Non-LTE methods
        NLTEMethod = inp_df[col0.isin(['NLTEMethod'])][1].values[0].upper()[0]
        if NLTEMethod == 'T':
            NLTEPath = 'None'
            vib_label = inp_df[col0.isin(['QNsVibLabel'])][1].iloc[0].replace(' ','').split(',')
            rot_label = inp_df[col0.isin(['QNsRotLabel'])][1].iloc[0].replace(' ','').split(',')
            Tvib_str = inp_df[col0.isin(['Tvib'])][1].values[0]
            Tvib_list = [int(tvib) for tvib in parse_TP_values(Tvib_str)]
            Trot_str = inp_df[col0.isin(['Trot'])][1].values[0]
            Trot_list = [int(trot) for trot in parse_TP_values(Trot_str)]
            # Align Tvib_list and Trot_list lengths for 2T NLTE:
            # - If one list has length 1, broadcast it to match the other.
            # - If both have length > 1, form the Cartesian product so that
            #   len(Tvib_list) == len(Trot_list) == len(original_Tvib) * len(original_Trot).
            if len(Tvib_list) == 1 and len(Trot_list) > 1:
                Tvib_list = Tvib_list * len(Trot_list)
            elif len(Tvib_list) > 1 and len(Trot_list) == 1:
                Trot_list = Trot_list * len(Tvib_list)
            elif len(Tvib_list) > 1 and len(Trot_list) > 1:
                Tvib_grid, Trot_grid = np.meshgrid(Tvib_list, Trot_list, indexing='ij')
                Tvib_list = Tvib_grid.ravel().tolist()
                Trot_list = Trot_grid.ravel().tolist()
            LTE_NLTE = '__2T-nlte'
        elif NLTEMethod == 'D':
            NLTEFile = inp_df[col0.isin(['NLTEMethod'])][2].values[0]
            NLTEPath = read_path + '/'.join(data_info) + '/' + NLTEFile
            NLTE_file_exist = os.path.exists(NLTEPath)   
            if NLTE_file_exist == False:
                NLTEPath = inp_df[col0.isin(['NLTEMethod'])][2].values[0]   
            else:
                pass   
            vib_label = []
            rot_label = []
            Tvib_list = []
            Trot_str = inp_df[col0.isin(['Trot'])][1].values[0]
            Trot_list = [int(trot) for trot in parse_TP_values(Trot_str)]
            # Align T_list length to match Trot_list for NLTE D
            if len(T_list) == 1 and len(Trot_list) > 1:
                T_list = T_list * len(Trot_list)
            elif len(T_list) > 1 and len(Trot_list) == 1:
                Trot_list = Trot_list * len(T_list)
            elif len(T_list) > 1 and len(Trot_list) > 1:
                T_grid, Trot_grid = np.meshgrid(T_list, Trot_list, indexing='ij')
                T_list = T_grid.ravel().tolist()
                Trot_list = Trot_grid.ravel().tolist()
            LTE_NLTE = '__nvib-nlte'
        elif NLTEMethod == 'P':
            NLTEFile = inp_df[col0.isin(['NLTEMethod'])][2].values[0]
            NLTEPath = read_path + '/'.join(data_info) + '/' + NLTEFile
            NLTE_file_exist = os.path.exists(NLTEPath)   
            if NLTE_file_exist == False:
                NLTEPath = inp_df[col0.isin(['NLTEMethod'])][2].values[0]   
            else:
                pass    
            vib_label = []
            rot_label = []
            Tvib_list = []
            Trot_list = []     
            LTE_NLTE = '__pop-nlte'
        else:
            raise ValueError("Please choose NLTE method 'TvibTrot' or 'Density' or 'Population'.")
    else:
        vib_label = []
        rot_label = []
        Tvib_list = []
        Trot_list = []
        NLTEMethod = 'L'
        NLTEPath = 'None'
        LTE_NLTE = ''

    # Stick spectra
    if StickSpectra != 0:
        PlotStickSpectraYN = inp_df[col0.isin(['PlotStickSpectra(Y/N)'])][1].values[0].upper()[0]
        PlotStickSpectraMethod = inp_df[col0.isin(['PlotStickSpectraMethod'])][1].values[0].upper()
        PlotStickSpectraWnWl = inp_df[col0.isin(['PlotStickSpectraWnWl'])][1].values[0].upper()
        PlotStickSpectraUnit = inp_df[col0.isin(['PlotStickSpectraWnWl'])][2].values[0].lower()
        if PlotStickSpectraYN == 'Y':
            _limitYaxisStickSpectra = inp_df[col0.isin(['Y-axisLimitStickSpectra'])][1].values[0]
            if _limitYaxisStickSpectra == '#':
                limitYaxisStickSpectra = 1e-30
            elif pd.isnull(_limitYaxisStickSpectra) == True:
                limitYaxisStickSpectra = 1e-30
            else:
                limitYaxisStickSpectra = float(_limitYaxisStickSpectra)
        else:
            limitYaxisStickSpectra = 0
    else:
        PlotStickSpectraYN = 'None'
        PlotStickSpectraMethod = 'None'
        PlotStickSpectraWnWl = 'None'
        PlotStickSpectraUnit = 'None'
        limitYaxisStickSpectra = 0
        
    # Cross sections
    if CrossSections != 0:
        NpointsORBinSize = inp_df[col0.isin(['Npoints/BinSize'])][1].values[0].upper()
        if 'POI' in NpointsORBinSize:
            N_point = int(inp_df[col0.isin(['Npoints/BinSize'])][2].iloc[0])
            bin_size = float(np.abs(max_wnl - min_wnl)/(N_point-1))
        elif 'BIN' in NpointsORBinSize or 'SIZ' in NpointsORBinSize:
            bin_size = float(inp_df[col0.isin(['Npoints/BinSize'])][2].iloc[0])
            N_point = int(np.abs(max_wnl - min_wnl)/bin_size+1)
        else:
            raise ValueError("Please type the correct grid choice 'Npoints' or 'BinSize' into the input file.")
        # Predissociative cross sections
        predissocYN = inp_df[col0.isin(['PredissocXsec(Y/N)'])][1].values[0].upper()[0]
        if predissocYN != 'Y' and predissocYN != 'N':
            raise ValueError("Please type the correct predissociative choice 'Y' or 'N' into the input file.")      
        if predissocYN == 'Y':
            photo = '__photo'
        else:
            photo = ''  
        # Cutoff
        cutoffYN = inp_df[col0.isin(['Cutoff(Y/N)'])][1].values[0].upper()[0]
        if cutoffYN == 'Y':
            cutoff = float(inp_df[col0.isin(['Cutoff(Y/N)'])][2].iloc[0])
        elif cutoffYN == 'N':
            cutoff = 'None'
        else:
            raise ValueError("Please type the correct cutoff choice 'Y' or 'N' into the input file.")
        # Other parameters
        broadeners = list(inp_df[col0.isin(['Broadeners'])].iloc[0])[1:]
        broadeners = [i for i in broadeners if i is not np.nan]
        ratios = np.array(list(inp_df[col0.isin(['Ratios'])].iloc[0])[1:], dtype=float)
        ratios = ratios[~np.isnan(ratios)]
        profile = inp_df[col0.isin(['Profile'])][1].values[0].upper().replace('PRO','')
        # Parse multiple pressure values
        if 'DOP' not in profile and 'GAU' not in profile:
            press_str = inp_df[col0.isin(['Pressures'])][1].values[0]
            P_list = [float(press) for press in parse_TP_values(press_str)]
        else:
            P_list = []
        # Doppler HWHM
        DopplerHWHMYN = inp_df[col0.isin(['DopplerHWHM(Y/N)'])][1].values[0].upper()[0]        
        if 'DOP' in profile: 
            alpha_HWHM = 'None'
            alpha_hwhm_colid = 0
        elif 'GAU' in profile:
            if DopplerHWHMYN == 'Y':
                alpha_HWHM = float(inp_df[col0.isin(['DopplerHWHM(Y/N)'])][2].iloc[0])
                alpha_hwhm_colid = 0
            elif DopplerHWHMYN == 'U':
                alpha_HWHM = []
                alpha_hwhm_colid = int(inp_df[col0.isin(['DopplerHWHM(Y/N)'])][2].iloc[0])
            else:
                raise ValueError("Gaussian line profile requires a HWHM. " 
                                  + "Please choose 'Y' and give a value for Doppler HWHM in the input file. " 
                                  + "Or 'U' and add your own Doppler HWHM values as the last column in the transitions file(s). "
                                  + "Otherwise, please choose Doppler line profile " 
                                  + "(with calculated temperature-dependent Doppler HWHM).")
        elif 'VOI' in profile:
            if DopplerHWHMYN == 'Y':
                alpha_HWHM = float(inp_df[col0.isin(['DopplerHWHM(Y/N)'])][2].iloc[0])
                alpha_hwhm_colid = 0
            elif DopplerHWHMYN == 'U':
                alpha_HWHM = []
                alpha_hwhm_colid = int(inp_df[col0.isin(['DopplerHWHM(Y/N)'])][2].iloc[0])
            elif DopplerHWHMYN == 'N':
                alpha_HWHM = 'None'
                alpha_hwhm_colid = 0
            else:
                raise ValueError("Please type the correct Doppler HWHM choice 'Y', 'N', or 'U' into the input file.")
        else:
            alpha_HWHM = 'None'
            alpha_hwhm_colid = 0
        # Lorentzian HWHM 
        LorentzianHWHMYN = inp_df[col0.isin(['LorentzianHWHM(Y/N)'])][1].values[0].upper()[0]  
        if LorentzianHWHMYN == 'Y':
            gamma_HWHM = float(inp_df[col0.isin(['LorentzianHWHM(Y/N)'])][2].iloc[0])
            gamma_hwhm_colid = 0
        elif LorentzianHWHMYN == 'U':
            gamma_HWHM = []
            gamma_hwhm_colid = int(inp_df[col0.isin(['LorentzianHWHM(Y/N)'])][2].iloc[0])
        elif LorentzianHWHMYN == 'N':
            gamma_HWHM = 'None'
            gamma_hwhm_colid = 0
        else:
            raise ValueError("Please type the correct Lorentzian HWHM choice 'Y', 'N', or 'U' into the input file.")
        # Plot 
        PlotCrossSectionYN = inp_df[col0.isin(['PlotCrossSection(Y/N)'])][1].values[0].upper()[0]  
        if PlotCrossSectionYN == 'Y':
            _limitYaxisXsec = inp_df[col0.isin(['Y-axisLimitXsec'])][1].values[0]
            if _limitYaxisXsec == '#':
                limitYaxisXsec = 1e-30
            elif pd.isnull(_limitYaxisXsec) == True:
                limitYaxisXsec = 1e-30
            else:
                limitYaxisXsec = float(_limitYaxisXsec)
        else:
            limitYaxisXsec = 0
        PlotCrossSectionMethod = inp_df[col0.isin(['PlotCrossSectionMethod'])][1].values[0].upper()
        PlotCrossSectionWnWl = inp_df[col0.isin(['PlotCrossSectionWnWl'])][1].values[0].upper()
        PlotCrossSectionUnit = inp_df[col0.isin(['PlotCrossSectionWnWl'])][2].values[0].lower() 
        if 'L' not in wn_wl:
            wn_grid = np.linspace(min_wnl, max_wnl, N_point)
        else:
            if wn_wl_unit == 'um':
                unit_change = 1e4
            elif wn_wl_unit == 'nm':
                unit_change = 1e7
            else:
                raise ValueError("Please wirte the unit of wavelength in the input file: um or nm.")
            if min_wnl == 0:
                raise ValueError("Please set the minimum wavenumber greater than 0 in the input file.")
            wl_grid = np.linspace(min_wnl, max_wnl, N_point)
            wn_grid = unit_change / wl_grid

    else:
        bin_size = 'None'
        N_point = 'None'
        predissocYN = 'N'
        photo = ''
        cutoff = 'None'   
        DopplerHWHMYN = 'None'
        LorentzianHWHMYN = 'None'      
        alpha_HWHM = 'None'        
        gamma_HWHM = 'None'
        alpha_hwhm_colid = 0
        gamma_hwhm_colid = 0
        broadeners = []
        ratios = np.array([])
        P_list = []
        wn_grid = np.linspace(0,1,1)
        profile = 'None'
        PlotCrossSectionYN = 'None'
        PlotCrossSectionMethod = 'None'
        PlotCrossSectionWnWl = 'None'
        PlotCrossSectionUnit = 'None'
        limitYaxisXsec = 0

    # Molecule and isotopologue ID, abundance, mass uncertainty, lifetime and g-factor           
    if database == 'ExoMol':
        # Read ExoMol definition file (.def.json).
        try:
            deffile_path = read_path + '/'.join(data_info) + '/' + '__'.join(data_info[-2:]) + '.def.json'  
            def_df = pd.read_json(deffile_path, orient='columns')
            states_file_fields = def_df['dataset']['states']['states_file_fields']
            states_file_fields_num = len(states_file_fields)
            states_col = [states_file_fields[i]['name'] for i in range(states_file_fields_num)]
            states_fmt = [states_file_fields[i]['cfmt'] for i in range(states_file_fields_num)]
            abundance = 1
            mass = def_df['isotopologue']['mass_in_Da']     # ExoMol mass (Dalton)
            try:
                # check_uncertainty = def_df['dataset']['states']['uncertainty_description']
                check_uncertainty = def_df['dataset']['states']['uncertainties_available']
            except:
                check_uncertainty = False
            try:
                check_predissoc = def_df['dataset']['predis'] or False
            except:
                check_predissoc = False
            try:
                check_lifetime = def_df['dataset']['states']['lifetime_available']
            except:
                check_lifetime = False
            try:
                check_gfactor = def_df['dataset']['states']['lande_g_available']
            except:
                check_gfactor = False   
        # Read ExoMol definition file (.def).
        except:
            deffile_path = read_path + '/'.join(data_info) + '/' + '__'.join(data_info[-2:]) + '.def'
            def_df = pd.read_csv(deffile_path,sep='\\s+',usecols=[0,1,2,3,4],names=['0','1','2','3','4'],header=None)
            states_cfmt_df = def_df[def_df['3'].isin(['Format'])]
            states_col = def_df.iloc[states_cfmt_df.index-1]['0'].tolist()
            states_fmt = states_cfmt_df['1'].tolist()
            abundance = 1
            mass = float(def_df[def_df['4'].isin(['mass'])]['0'].values[0])     # ExoMol mass (Dalton)
            try:
                check_uncertainty = bool(int(def_df[def_df['2'].isin(['Uncertainty'])]['0'].values[0]))
            except:
                check_uncertainty = False
            try:
                check_predissoc = bool(int(def_df[def_df['2'].isin(['Predissociative'])]['0'].values[0]))
            except:
                check_predissoc = False
            try:
                check_lifetime = bool(int(def_df[def_df['2'].isin(['Lifetime'])]['0'].values[0]))
            except:
                check_lifetime = False
            try:
                check_gfactor = bool(int(def_df[def_df['3'].isin(['g-factor'])]['0'].values[0]))
            except:
                check_gfactor = False

    elif database == 'ExoAtom':
        # Read ExoAtom definition file (.adef.json).
        deffile_path = read_path + '/'.join(data_info) + '/' + '__'.join(data_info[-2:]) + '.adef.json'  
        def_df = pd.read_json(deffile_path, orient='columns')
        states_file_fields = def_df['dataset']['states']['states_file_fields']
        states_file_fields_num = len(states_file_fields)
        states_col = [states_file_fields[i]['name'] for i in range(states_file_fields_num)]
        states_fmt = [states_file_fields[i]['cfmt'] for i in range(states_file_fields_num)]
        abundance = 1
        mass = def_df['species']['mass_in_Da']     # ExoAtom mass (Dalton)
        try:
            # check_uncertainty = def_df['dataset']['states']['uncertainty_description']
            check_uncertainty = def_df['dataset']['states']['uncertainty_available']
        except:
            check_uncertainty = False
        try:
            check_predissoc = def_df['dataset']['predis'] or False
        except:
            check_predissoc = False
        try:
            check_lifetime = def_df['dataset']['states']['lifetime_available']
        except:
            check_lifetime = False
        try:
            check_gfactor = def_df['dataset']['states']['lande_g_available']
        except:
            check_gfactor = False

    elif database == 'HITRAN' or database == 'HITEMP':
        read_path = read_path.rstrip('/')
        try:
            # read_path = /.../Databases/HITRAN/NO/NO__14N-16O.par
            # molparam_filepath = /.../Databases/HITRAN/molparam.txt
            hitran_root = os.path.dirname(os.path.dirname(read_path))
            molparam_filepath = os.path.join(hitran_root, 'molparam.txt')
            # Parse molparam.txt to extract abundance and molar mass for this species_main_id
            # and construct a DataFrame with the same columns used below:
            # 'local ID', 'Abundance', 'Molar Mass /g·mol-1'
            local_ids = []
            abundances = []
            molar_masses = []
            current_mol_id = None
            try:
                with open(molparam_filepath, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        # Skip header line
                        if line.startswith('Molecule #'):
                            continue
                        # Molecule header lines look like: "   NO (8)"
                        m = re.search(r'\((\d+)\)\s*$', line)
                        if m:
                            current_mol_id = int(m.group(1))
                            continue
                        # Only process isotopologue lines for the target molecule
                        if current_mol_id != species_main_id:
                            continue
                        # Data lines for isotopologues start with an integer ISO code
                        parts = line.split()
                        if len(parts) < 5:
                            continue
                        try:
                            # parts[1] = Abundance, parts[4] = Molar Mass(g)
                            abundance_val = float(parts[1])
                            molar_mass_val = float(parts[4])
                        except ValueError:
                            continue
                        local_ids.append(len(local_ids) + 1)  # 1-based local ID ordering
                        abundances.append(abundance_val)
                        molar_masses.append(molar_mass_val)
            except FileNotFoundError:
                raise FileNotFoundError(
                    f"Failed to load HITRAN iso-meta data: neither online table nor local molparam.txt found. "
                    f"Expected molparam at: {molparam_filepath}"
                )

            if not local_ids:
                raise RuntimeError(
                    f"No isotopologue entries found in molparam.txt for HITRAN molecule ID {species_main_id} "
                    f"(molparam file: {molparam_filepath})"
                )

            iso_meta_table = pd.DataFrame({
                'local ID': np.array(local_ids, dtype=int),
                'Abundance': np.array(abundances, dtype=float),
                'Molar Mass /g·mol-1': np.array(molar_masses, dtype=float),
            })
        except Exception:
            # Primary: try to read iso-meta table from HITRAN website
            isometa_url = 'https://hitran.org/docs/iso-meta/'
            iso_meta_table = pd.read_html(isometa_url)[species_main_id - 1]

        iso_meta_row = iso_meta_table[iso_meta_table['local ID'].isin([species_sub_id])]
        states_col = []
        states_fmt = []
        # Online table stores Abundance as string with special formatting, while
        # the local molparam-based table stores it as a plain float. Handle both.
        abundance_val = iso_meta_row['Abundance'].iloc[0]
        if isinstance(abundance_val, str):
            abundance = float(abundance_val.replace('\xa0×\xa010','E'))
        else:
            abundance = float(abundance_val)
        mass = float(iso_meta_row['Molar Mass /g·mol-1'].iloc[0])                   # HITRAN molar mass (g/mol)
        check_uncertainty = True
        check_predissoc = False
        check_lifetime = False
        check_gfactor = False

    else:
        raise ValueError("Please add the name of the database 'ExoMol', 'ExoAtom', 'HITRAN', or 'HITEMP' into the input file.")

    if predissocYN == 'Y' and check_predissoc is False and 'LOR' not in profile:
        print(
            "Predissociation cross sections are enabled (PredissocXsec='Y') and a non-Lorentzian profile is selected.\n"
            "The program may take a long time because predissociative lifetimes must be calculated before the cross sections.\n"
        )
            
    return (database, data_info, read_path, save_path, logs_path,
            Conversion, PartitionFunctions, SpecificHeats, CoolingFunctions, Lifetimes, OscillatorStrengths, StickSpectra, CrossSections,
            ncputrans, ncpufiles, chunk_size, ConversionFormat, ConversionMinFreq, ConversionMaxFreq, ConversionUnc, ConversionThreshold, 
            GlobalQNLabel_list, GlobalQNFormat_list, LocalQNLabel_list, LocalQNFormat_list,
            Ntemp, Tmax, CompressYN, gfORf, broadeners, ratios, T_list, P_list, wn_wl, wn_wl_unit, min_wnl, max_wnl, min_wn, max_wn, N_point, bin_size, wn_grid, 
            predissocYN, photo, cutoff, threshold, UncFilter, NLTEMethod, NLTEPath, LTE_NLTE, QNslabel_list, QNsformat_list, QNs_label, QNs_value, QNs_format, QNsFilter, 
            DopplerHWHMYN, LorentzianHWHMYN, alpha_HWHM, gamma_HWHM, alpha_hwhm_colid, gamma_hwhm_colid, abs_emi, profile, 
            species_id, species_main_id, species_sub_id, abundance, mass, states_col, states_fmt,
            check_uncertainty, check_lifetime, check_gfactor, check_predissoc, 
            PlotOscillatorStrengthYN, PlotOscillatorStrengthMethod, PlotOscillatorStrengthWnWl, PlotOscillatorStrengthUnit, limitYaxisOS,
            PlotStickSpectraYN, PlotStickSpectraMethod, PlotStickSpectraWnWl, PlotStickSpectraUnit, limitYaxisStickSpectra, 
            Tvib_list, Trot_list, vib_label, rot_label, PlotCrossSectionYN, PlotCrossSectionMethod, PlotCrossSectionWnWl, PlotCrossSectionUnit, limitYaxisXsec)

