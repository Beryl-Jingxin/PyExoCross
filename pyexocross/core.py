"""
Core execution module for PyExoCross.

Contains the main get_results function that orchestrates all calculations.
"""
import os
import sys
import glob
import numpy as np
from tabulate import tabulate
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

# Add src directory to path for imports
# We need to add the parent directory (project root) to sys.path so that 'src' is recognized as a package
_current_dir = os.path.dirname(os.path.abspath(__file__))
# _current_dir is pyexocross/, so parent is project root
_project_root = os.path.dirname(_current_dir)
if _project_root not in sys.path and os.path.exists(_project_root):
    sys.path.insert(0, _project_root)

# Import all necessary modules using absolute imports from src package
# Note: We import at function level to avoid circular import issues
# These will be imported when needed in get_results()


def get_results(config):
    """
    Main function to execute PyExoCross calculations.

    Parameters
    ----------
    config : Config
        Configuration object containing all parameters.
    """
    # Import modules here to avoid circular import issues
    from src.base.utils import Timer
    from src.database.load_exomol import read_all_states, read_part_states, get_transfiles
    from src.database.load_hitran import read_parfile, read_hitran_parfile
    from src.process.Q_multi_T import cal_pf_multiT
    from src.process.hitran_qn import hitran_linelist_QN
    from src.calculation.calculate_lifetime import cal_lifetime
    from src.base.large_file import read_trans_chunks
    from src.conversion.exomol_to_hitran import conversion_exomol2hitran
    from src.conversion.hitran_to_exomol import conversion_hitran2exomol
    from src.save.exomol.exomol_partition_func import save_exomol_partition_func
    from src.save.exomol.exomol_specific_heat import save_exomol_specific_heat
    from src.save.exomol.exomol_lifetime import save_exomol_lifetime
    from src.save.exomol.exomol_cooling_func import save_exomol_cooling_func
    from src.save.exomol.exomol_oscillator_strength import save_exomol_oscillator_strength
    from src.save.exomol.exomol_stick_spectra import save_exomol_stick_spectra
    from src.save.exomol.exomol_cross_section import save_exomol_cross_section
    from src.save.exomol.exomol_stick_spectra_cross_section import save_exomol_stick_spectra_cross_section
    from src.save.hitran.hitran_partition_func import save_hitran_partition_func
    from src.save.hitran.hitran_specific_heat import save_hitran_specific_heat
    from src.save.hitran.hitran_lifetime import save_hitran_lifetime
    from src.save.hitran.hitran_cooling_func import save_hitran_cooling_func
    from src.save.hitran.hitran_oscillator_strength import save_hitran_oscillator_strength
    from src.save.hitran.hitran_stick_spectra import save_hitran_stick_spectra
    from src.save.hitran.hitran_cross_section import save_hitran_cross_section
    from src.save.hitran.hitran_stick_spectra_cross_section import save_hitran_stick_spectra_cross_section
    
    # Set globals from config (for legacy-style modules)
    config.to_globals()
    
    # Update bin-size–dependent constants used by line profile and
    # cross-section routines so that modules importing from
    # src.base.constants see consistent values.
    from src.base import constants as _const_mod
    bin_const_dict = _const_mod.get_bin_size_constants(config.bin_size)
    for _name, _val in bin_const_dict.items():
        setattr(_const_mod, _name, _val)
    
    # Access globals (set by config.to_globals())
    database = config.database
    data_info = config.data_info
    read_path = config.read_path
    save_path = config.save_path
    
    # Function flags
    Conversion = config.conversion
    PartitionFunctions = config.partition_functions
    SpecificHeats = config.specific_heats
    Lifetimes = config.lifetimes
    CoolingFunctions = config.cooling_functions
    OscillatorStrengths = config.oscillator_strengths
    StickSpectra = config.stick_spectra
    CrossSections = config.cross_sections
    
    # Other parameters
    Ntemp = config.ntemp
    Tmax = config.tmax
    T_list = config.T_list
    P_list = config.P_list
    Tvib_list = config.tvib_list
    Trot_list = config.trot_list
    ConversionFormat = config.conversion_format
    ConversionMinFreq = config.conversion_min_freq
    ConversionMaxFreq = config.conversion_max_freq
    ConversionUnc = config.conversion_unc
    ConversionThreshold = config.conversion_threshold
    min_wn = config.min_wn
    max_wn = config.max_wn
    UncFilter = config.unc_filter
    threshold = config.threshold
    profile = config.profile
    predissocYN = config.predissoc_yn
    check_predissoc = config.check_predissoc
    check_lifetime = config.check_lifetime
    states_col = config.states_col
    states_fmt = config.states_fmt
    ncpufiles = config.ncpufiles
    ncputrans = config.ncputrans
    
    t_tot = Timer()
    t_tot.start()
    
    # Print database information
    if database == 'ExoMol':
        print('ExoMol database')
        headers = ['Molecule', 'Isotopologue', 'Dataset']
        print(tabulate([data_info], headers=headers, tablefmt="fancy_grid"))
    elif database == 'ExoAtom':
        print('ExoAtom database')
        headers = ['Atom', 'Dataset']
        print(tabulate([data_info], headers=headers, tablefmt="fancy_grid"))
    elif database == 'HITRAN' or database == 'HITEMP':
        print(database, 'database')
        headers = ['Molecule', 'Molecule ID', 'Isotopologue', 'Isotopologue ID', 'Dataset']
        species_main_id = config.species_main_id
        species_sub_id = config.species_sub_id
        print(tabulate([[data_info[0], str(species_main_id), data_info[1], str(species_sub_id), data_info[2]]],
                       headers=headers, tablefmt="fancy_grid"))
    
    if database == 'ExoMol' or database == 'ExoAtom':
        # Functions
        Nfunctions = (Conversion + PartitionFunctions + SpecificHeats + Lifetimes
                     + CoolingFunctions + OscillatorStrengths + StickSpectra + CrossSections)
        if Nfunctions > 0:
            states_df = read_all_states(
                read_path,
                data_info,
                config.check_uncertainty,
                config.states_col,
                config.states_fmt,
            )
        else:
            raise ValueError("Please choose functions which you want to calculate.")
        
        # These functions need whole states
        NeedAllStates = (Conversion + PartitionFunctions + SpecificHeats
                        + Lifetimes + CoolingFunctions + OscillatorStrengths)
        if NeedAllStates > 0:
            if (Conversion == 1 and ConversionFormat == 1):
                conversion_exomol2hitran(states_df)
            if PartitionFunctions == 1:
                save_exomol_partition_func(states_df, Ntemp, Tmax)
            if SpecificHeats == 1:
                save_exomol_specific_heat(states_df, Ntemp, Tmax)
            if Lifetimes == 1:
                save_exomol_lifetime(read_path, states_df, states_col, states_fmt)
            if CoolingFunctions == 1:
                save_exomol_cooling_func(states_df, Ntemp, Tmax)
            if OscillatorStrengths == 1:
                save_exomol_oscillator_strength(states_df)
        
        # Only calculating stick spectra and cross sections need part of states
        NeedPartStates = StickSpectra + CrossSections
        if NeedPartStates > 0:
            states_part_df = read_part_states(
                states_df,
                config.unc_filter,
                config.nlte_method,
                config.nlte_path,
                config.check_uncertainty,
                config.check_lifetime,
                config.check_gfactor,
                config.qnslabel_list,
                config.qns_label,
                config.qns_filter,
                config.qns_value,
                config.vib_label,
                config.rot_label,
            )
            Q_arr = cal_pf_multiT(
                T_list,
                Tvib_list,
                Trot_list,
                states_part_df,
                config.nlte_method,
                config.read_path,
                config.data_info,
            )
            states_df.index.name = 'index'
            # Calculate predissociation spectra and cross sections if lifetimes are not exist in the states file
            if predissocYN == 'Y' and check_predissoc + check_lifetime == 0 and 'VOI' in profile:
                np.seterr(divide='ignore', invalid='ignore')
                print('Calculate lifetimes.')
                t = Timer()
                t.start()
                print('Reading all transitions and calculating lifetimes ...')
                trans_filepaths = get_transfiles(read_path, data_info)
                # For each transitions file, stream chunks and run cal_lifetime in parallel.
                # Keep cal_lifetime API unchanged: it takes (states_df, trans_df_chunk).
                # Use the full states_df for lifetime accumulation so transition indices
                # from chunks always fall within valid state-id bounds; then map the
                # resulting tau values back to states_part_df.
                lifetime_sum = np.zeros(len(states_df), dtype=float)
                for trans_filepath in trans_filepaths:
                    trans_filename = trans_filepath.split('/')[-1]
                    print('Processing transitions file for predissociation lifetime:', trans_filename)
                    use_cols = [0, 1, 2]
                    use_names = ['u', 'l', 'A']
                    trans_reader = read_trans_chunks(trans_filepath, use_cols, use_names)
                    trans_chunks = list(trans_reader)
                    if len(trans_chunks) == 0:
                        continue
                    with ProcessPoolExecutor(max_workers=ncputrans) as trans_executor:
                        futures = [
                            trans_executor.submit(cal_lifetime, states_df, trans_chunk)
                            for trans_chunk in tqdm(trans_chunks, desc='Processing ' + trans_filename)
                        ]
                        lifetime_sum += np.sum(
                            [future.result() for future in tqdm(futures, desc='Combining ' + trans_filename)],
                            axis=0,
                        )
                # Avoid division by zero; entries with zero A-sum correspond to infinite lifetime
                with np.errstate(divide='ignore', invalid='ignore'):
                    lifetime_result = np.where(lifetime_sum > 0, 1.0 / lifetime_sum, np.inf)
                states_df_tau = np.array(
                    [f'{x:>12.4E}'.replace('INF', 'Inf') for x in lifetime_result]
                ).astype(float)
                # Map full-states tau values onto the filtered states_part_df by state id.
                tau_map = dict(zip(states_df['id'].values, states_df_tau))
                states_part_df['tau'] = states_part_df['id'].map(tau_map).astype(float)
                t.end()
                print('Finished reading all transitions and calculating lifetimes!\n')
            # If both stick spectra and cross sections are needed, process together to read files once
            if StickSpectra == 1 and CrossSections == 1:
                save_exomol_stick_spectra_cross_section(
                    states_part_df,
                    T_list,
                    Tvib_list,
                    Trot_list,
                    P_list,
                    Q_arr,
                    config.check_uncertainty,
                    config.check_lifetime,
                    config.check_gfactor,
                )
            elif StickSpectra == 1:
                save_exomol_stick_spectra(states_part_df, T_list, Tvib_list, Trot_list, Q_arr)
            elif CrossSections == 1:
                save_exomol_cross_section(states_part_df, T_list, Tvib_list, Trot_list, P_list, Q_arr)
    
    elif database == 'HITRAN' or database == 'HITEMP':
        parfile_df = read_parfile(read_path)
        if (Conversion == 1 and ConversionFormat == 2):
            hitran_df = read_hitran_parfile(read_path,parfile_df,ConversionMinFreq,ConversionMaxFreq,
                                            ConversionUnc,ConversionThreshold).reset_index().drop(columns='index')
            (hitran_states_col, hitran_states_fmt) = conversion_hitran2exomol(hitran_df)
        # Use ExoMol functions
        NuseExoMolFunc = PartitionFunctions+SpecificHeats+Lifetimes
        if NuseExoMolFunc > 0:
            read_hitran2exomol_path = save_path + 'conversion/HITRAN2ExoMol/'
            conversion_foldername = read_hitran2exomol_path+'/'.join(data_info)+'/'
            if os.path.exists(conversion_foldername):
                states_list = glob.glob(conversion_foldername + '__'.join(data_info[-2:]) + '.states.bz2')
                trans_list = glob.glob(conversion_foldername + '__'.join(data_info[-2:]) + '*.trans.bz2')
                hitran_df = read_hitran_parfile(read_path,parfile_df,ConversionMinFreq,ConversionMaxFreq,
                                                'None','None').reset_index().drop(columns='index')
                (hitran_states_col, hitran_states_fmt) = conversion_hitran2exomol(hitran_df)
            # These functions need whole states
            states_df = read_all_states(
                read_hitran2exomol_path,
                data_info,
                config.check_uncertainty,
                config.states_col,
                config.states_fmt,
            )
            if PartitionFunctions == 1:
                save_hitran_partition_func(states_df, Ntemp, Tmax)
            if SpecificHeats == 1:
                save_hitran_specific_heat(states_df, Ntemp, Tmax)
            if Lifetimes == 1:
                save_hitran_lifetime(read_hitran2exomol_path, states_df, hitran_states_col, hitran_states_fmt)
        # Use HITRAN functions
        # Calculate cooling functions or oscillator strengths
        NeedWholeHITRAN = CoolingFunctions + OscillatorStrengths
        if NeedWholeHITRAN != 0:
            hitran_df = read_hitran_parfile(read_path, parfile_df, min_wn, max_wn, 'None', 'None').reset_index().drop(columns='index')
        if CoolingFunctions == 1:
            save_hitran_cooling_func(hitran_df, Ntemp, Tmax)
        if OscillatorStrengths == 1:
            save_hitran_oscillator_strength(hitran_df)
        if (StickSpectra + CrossSections > 0):
            hitran_df = read_hitran_parfile(read_path, parfile_df, min_wn, max_wn,
                                           UncFilter, threshold).reset_index().drop(columns='index')
            (hitran_linelist_df, QNs_col) = hitran_linelist_QN(hitran_df)
        if StickSpectra == 1 and CrossSections == 1:
            # Use combined function for both stick spectra and cross sections
            save_hitran_stick_spectra_cross_section(hitran_linelist_df, QNs_col, T_list, P_list, Tvib_list, Trot_list)
        elif StickSpectra == 1:
            save_hitran_stick_spectra(hitran_linelist_df, QNs_col, T_list, Tvib_list, Trot_list)
        elif CrossSections == 1:
            save_hitran_cross_section(hitran_linelist_df, T_list, P_list, Tvib_list, Trot_list)
    else:
        raise ValueError("Please add the name of the database 'ExoMol', 'ExoAtom', 'HITRAN', or 'HITEMP' into the input file.")
    
    print('The program total running time:')
    t_tot.end()
    print('\nFinished!')

