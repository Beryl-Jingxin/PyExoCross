"""
Core execution module for PyExoCross.

Contains the main get_results function that orchestrates all calculations.
"""
import copy
import glob
import os
import numpy as np
from tabulate import tabulate
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor


def printdeviceinfo(config):
    """Print the selected CPU/GPU runtime table."""
    from pyexocross.gpu.base_gpu import configure_runtime

    configure_runtime(
        run_mode=getattr(config, 'run_mode', 'CPU'),
        gpu_backend=getattr(config, 'gpu_backend', 'AUTO'),
        gpu_batch_lines=getattr(config, 'gpu_batch_lines', 8192),
        gpu_batch_grid=getattr(config, 'gpu_batch_grid', 256),
        verbose=True,
    )


def printdatabaseinfo(config):
    """Print database and species metadata in one table."""
    database = config.database
    data_info = config.data_info
    if database in ('ExoMol', 'ExoMolHR'):
        headers = ['Database', 'Molecule', 'Isotopologue', 'Dataset']
        row = [database] + data_info
    elif database == 'ExoAtom':
        headers = ['Database', 'Atom', 'Dataset']
        row = [database] + data_info
    elif database in ('HITRAN', 'HITEMP'):
        headers = [
            'Database',
            'Molecule',
            'Molecule ID',
            'Isotopologue',
            'Isotopologue ID',
            'Dataset',
        ]
        row = [
            database,
            data_info[0],
            str(config.species_main_id),
            data_info[1],
            str(config.species_sub_id),
            data_info[2],
        ]
    else:
        return
    print(tabulate([row], headers=headers, tablefmt='fancy_grid'))


def get_results(config, data=None):
    """
    Main function to execute PyExoCross calculations.

    Parameters
    ----------
    config : Config
        Configuration object containing all parameters.
    """
    # Print run identity before importing calculation modules or starting timers.
    config.to_globals()
    if data is None:
        printdeviceinfo(config)
        printdatabaseinfo(config)

    # Import modules here to avoid circular import issues
    from pyexocross.base.utils import Timer
    from pyexocross.database.load_exomol import read_all_states, read_part_states, get_transfiles
    from pyexocross.database.load_exomolhr import read_exomolhr_df
    from pyexocross.database.load_hitran import read_parfile, read_hitran_parfile
    from pyexocross.process.Q_multi_T import cal_pf_multiT
    from pyexocross.process.filter_qn import QNfilter_linelist
    from pyexocross.process.hitran_qn import hitran_linelist_QN
    from pyexocross.calculation.calculate_lifetime import cal_lifetime
    from pyexocross.base.large_file import read_trans_chunks, sourcename
    from pyexocross.convert.exomol_to_hitran import conversion_exomol2hitran
    from pyexocross.convert.hitran_to_exomol import conversion_hitran2exomol
    from pyexocross.convert.exomolhr_to_hitran import conversion_exomolhr2hitran
    from pyexocross.save.exomol.exomol_partition_func import save_exomol_partition_func
    from pyexocross.save.exomol.exomol_specific_heat import save_exomol_specific_heat
    from pyexocross.save.exomol.exomol_lifetime import save_exomol_lifetime
    from pyexocross.save.exomol.exomol_cooling_func import save_exomol_cooling_func
    from pyexocross.save.exomol.exomol_oscillator_strength import save_exomol_oscillator_strength
    from pyexocross.save.exomol.exomol_stick_spectra import save_exomol_stick_spectra
    from pyexocross.save.exomol.exomol_cross_section import save_exomol_cross_section
    from pyexocross.save.exomol.exomol_stick_spectra_cross_section import save_exomol_stick_spectra_cross_section
    from pyexocross.save.exomolhr.exomolhr_stick_spectra import save_exomolhr_stick_spectra
    from pyexocross.save.exomolhr.exomolhr_cross_section import save_exomolhr_cross_section
    from pyexocross.save.exomolhr.exomolhr_stick_spectra_cross_section import save_exomolhr_stick_spectra_cross_section
    from pyexocross.save.hitran.hitran_partition_func import save_hitran_partition_func
    from pyexocross.save.hitran.hitran_specific_heat import save_hitran_specific_heat
    from pyexocross.save.hitran.hitran_lifetime import save_hitran_lifetime
    from pyexocross.save.hitran.hitran_cooling_func import save_hitran_cooling_func
    from pyexocross.save.hitran.hitran_oscillator_strength import save_hitran_oscillator_strength
    from pyexocross.save.hitran.hitran_stick_spectra import save_hitran_stick_spectra
    from pyexocross.save.hitran.hitran_cross_section import save_hitran_cross_section
    from pyexocross.save.hitran.hitran_stick_spectra_cross_section import save_hitran_stick_spectra_cross_section
    
    # Update bin-size–dependent constants used by line profile and
    # cross-section routines so that modules importing from
    # pyexocross.base.constants see consistent values.
    from pyexocross.base import constants as _const_mod
    bin_const_dict = _const_mod.get_bin_size_constants(config.bin_size)
    for _name, _val in bin_const_dict.items():
        setattr(_const_mod, _name, _val)
        
    # Pre-calculate Doppler HWHM constant outside the def and inject it
    # as a global constant for Doppler_HWHM to use directly.
    import pyexocross.calculation.calcualte_line_profile as _calc_mod
    _calc_mod.Sqrt2NAkBln2mInvc = _const_mod.get_doppler_constants(config.mass)
    
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
    
    if database == 'ExoMol' or database == 'ExoAtom':
        # Functions
        Nfunctions = (Conversion + PartitionFunctions + SpecificHeats + Lifetimes
                     + CoolingFunctions + OscillatorStrengths + StickSpectra + CrossSections)
        NeedPartStates = StickSpectra + CrossSections
        NeedAllStates = (Conversion + PartitionFunctions + SpecificHeats
                        + Lifetimes + CoolingFunctions + OscillatorStrengths)
        NeedConversionTransitions = int(Conversion == 1 and ConversionFormat == 'HITRAN')
        NeedAllTransitions = Lifetimes + CoolingFunctions + OscillatorStrengths
        if (
            predissocYN == 'Y'
            and check_predissoc + check_lifetime == 0
            and 'VOI' in profile
        ):
            NeedAllTransitions += 1
        UsesTransitions = NeedAllTransitions + NeedPartStates + NeedConversionTransitions
        if data is not None and NeedAllTransitions > 0 and not data.alltrans:
            raise ValueError(
                'This calculation needs every transition file. '
                'Call px.load(..., all_transitions=True).'
            )
        if data is None and UsesTransitions > 0 and config.cache != 'none':
            from pyexocross.database.data import loaddata

            loadconfig = config
            if NeedAllTransitions == 0 and NeedConversionTransitions > 0:
                loadconfig = copy.copy(config)
                if NeedPartStates > 0:
                    loadconfig.min_wn = min(config.min_wn, ConversionMinFreq)
                    loadconfig.max_wn = max(config.max_wn, ConversionMaxFreq)
                else:
                    loadconfig.min_wn = ConversionMinFreq
                    loadconfig.max_wn = ConversionMaxFreq
                loadconfig.cutoff = 'None'
            data = loaddata(
                loadconfig,
                cache=config.cache,
                cachedir=config.cache_dir,
                maxmemory=config.max_memory,
                refresh=config.refresh_cache,
                alltrans=NeedAllTransitions > 0,
                preparestates=NeedPartStates > 0,
            )
        if Nfunctions > 0:
            if data is None:
                states_df = read_all_states(
                    read_path,
                    data_info,
                    config.check_uncertainty,
                    config.states_col,
                    config.states_fmt,
                )
            elif NeedAllStates > 0:
                states_df = data.fullstates().copy()
            else:
                states_df = data.states.copy() if data.states is not None else None
        else:
            raise ValueError("Please choose functions which you want to calculate.")
        
        # These functions need whole states
        all_trans_sources = (
            data.transitions
            if data is not None and data.alltrans
            else None
        )
        if NeedAllStates > 0:
            if (Conversion == 1 and ConversionFormat == 'HITRAN'):
                conversion_sources = (
                    data.sources(ConversionMinFreq, ConversionMaxFreq)
                    if data is not None
                    else None
                )
                conversion_exomol2hitran(
                    states_df,
                    trans_sources=conversion_sources,
                )
            if PartitionFunctions == 1:
                save_exomol_partition_func(states_df, Ntemp, Tmax)
            if SpecificHeats == 1:
                save_exomol_specific_heat(states_df, Ntemp, Tmax)
            if Lifetimes == 1:
                save_exomol_lifetime(
                    read_path,
                    states_df,
                    states_col,
                    states_fmt,
                    trans_sources=all_trans_sources,
                )
            if CoolingFunctions == 1:
                save_exomol_cooling_func(
                    states_df,
                    Ntemp,
                    Tmax,
                    trans_sources=all_trans_sources,
                )
            if OscillatorStrengths == 1:
                save_exomol_oscillator_strength(
                    states_df,
                    trans_sources=all_trans_sources,
                )
        
        # Only calculating stick spectra and cross sections need part of states
        if NeedPartStates > 0:
            if data is None:
                states_part_df = read_part_states(
                    states_df,
                    config.unc_filter,
                    config.nlte_method,
                    config.nlte_path,
                    config.check_uncertainty,
                    config.check_lifetime,
                    config.check_gfactor,
                    config.qnslabel_list,
                    config.states_col,
                    config.states_fmt,
                    config.qns_label,
                    config.qns_filter,
                    config.qns_value,
                    config.vib_label,
                    config.rot_label,
                )
                trans_sources = None
            else:
                data.ensureqns(config)
                states_part_df = data.prepared(config)
                trans_sources = data.sources(
                    config.min_wn,
                    config.max_wn,
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
            if states_df is not None:
                states_df.index.name = 'index'
            # Calculate predissociation spectra and cross sections if lifetimes are not exist in the states file
            if predissocYN == 'Y' and check_predissoc + check_lifetime == 0 and 'VOI' in profile:
                np.seterr(divide='ignore', invalid='ignore')
                print('Calculate lifetimes.')
                t = Timer()
                t.start()
                print('Reading all transitions and calculating lifetimes ...')
                trans_filepaths = (
                    all_trans_sources
                    if all_trans_sources is not None
                    else get_transfiles(read_path, data_info)
                )
                # For each transitions file, stream chunks and run cal_lifetime in parallel.
                # Keep cal_lifetime API unchanged: it takes (states_df, trans_df_chunk).
                # Use the full states_df for lifetime accumulation so transition indices
                # from chunks always fall within valid state-id bounds; then map the
                # resulting tau values back to states_part_df.
                lifetime_sum = np.zeros(len(states_df), dtype=float)
                for trans_filepath in trans_filepaths:
                    trans_filename = sourcename(trans_filepath)
                    print('Processing transitions file for predissociation lifetime:', trans_filename)
                    use_cols = [0, 1, 2]
                    use_names = ['uid', 'lid', 'A']
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
                    trans_sources=trans_sources,
                )
            elif StickSpectra == 1:
                save_exomol_stick_spectra(
                    states_part_df,
                    T_list,
                    Tvib_list,
                    Trot_list,
                    Q_arr,
                    trans_sources=trans_sources,
                )
            elif CrossSections == 1:
                save_exomol_cross_section(
                    states_part_df,
                    T_list,
                    Tvib_list,
                    Trot_list,
                    P_list,
                    Q_arr,
                    trans_sources=trans_sources,
                )
    
    elif database == 'ExoMolHR':
        unsupported_funcs = PartitionFunctions + SpecificHeats + Lifetimes + CoolingFunctions + OscillatorStrengths
        if unsupported_funcs > 0:
            raise ValueError(
                "ExoMolHR currently supports conversion, stick spectra and cross sections only."
            )
        if predissocYN == 'Y':
            raise ValueError("Predissociation cross sections are not supported for ExoMolHR.")
        if Conversion + StickSpectra + CrossSections == 0:
            raise ValueError("Please choose conversion, stick spectra, or cross sections for ExoMolHR.")
        if data is None and config.cache != 'none':
            from pyexocross.database.data import loaddata

            data = loaddata(
                config,
                cache=config.cache,
                cachedir=config.cache_dir,
                maxmemory=config.max_memory,
                refresh=config.refresh_cache,
            )
        
        # Conversion: ExoMolHR → HITRAN
        if Conversion == 1 and ConversionFormat == 'HITRAN':
            if data is None:
                conv_exomolhr_df = read_exomolhr_df(
                    read_path,
                    data_info,
                    config.conversion_min_freq,
                    config.conversion_max_freq,
                    config.conversion_unc,
                )
            else:
                conv_exomolhr_df = data.lines(
                    config.conversion_min_freq,
                    config.conversion_max_freq,
                    config.conversion_unc,
                )
            conversion_exomolhr2hitran(conv_exomolhr_df)
        
        if StickSpectra + CrossSections > 0:
            if data is None:
                exomolhr_df = read_exomolhr_df(
                    read_path,
                    data_info,
                    min_wn,
                    max_wn,
                    UncFilter,
                )
            else:
                exomolhr_df = data.lines(min_wn, max_wn, UncFilter)
            if config.qns_filter != []:
                exomolhr_df = QNfilter_linelist(exomolhr_df, config.qns_value, config.qns_label)
            QNs_col = [label + "'" for label in config.qnslabel_list] + [label + '"' for label in config.qnslabel_list]

            if StickSpectra == 1 and CrossSections == 1:
                save_exomolhr_stick_spectra_cross_section(
                exomolhr_df,
                QNs_col,
                T_list,
                Tvib_list,
                Trot_list,
                P_list,
                )
            elif StickSpectra == 1:
                save_exomolhr_stick_spectra(
                    exomolhr_df,
                    QNs_col,
                    T_list,
                    Tvib_list,
                    Trot_list,
                )
            elif CrossSections == 1:
                save_exomolhr_cross_section(
                    exomolhr_df,
                    T_list,
                    P_list,
                    Tvib_list,
                    Trot_list,
                )

    elif database == 'HITRAN' or database == 'HITEMP':
        if data is None and config.cache != 'none':
            from pyexocross.database.data import loaddata

            data = loaddata(
                config,
                cache=config.cache,
                cachedir=config.cache_dir,
                maxmemory=config.max_memory,
                refresh=config.refresh_cache,
            )
        parfile_df = read_parfile(read_path) if data is None else None
        if (Conversion == 1 and ConversionFormat == 'ExoMol'):
            if data is None:
                hitran_df = read_hitran_parfile(
                    read_path,
                    parfile_df,
                    ConversionMinFreq,
                    ConversionMaxFreq,
                    ConversionUnc,
                    ConversionThreshold,
                ).reset_index(drop=True)
            else:
                hitran_df = data.lines(
                    ConversionMinFreq,
                    ConversionMaxFreq,
                    ConversionUnc,
                    ConversionThreshold,
                )
            (hitran_states_col, hitran_states_fmt) = conversion_hitran2exomol(hitran_df)
        # Use ExoMol functions
        NuseExoMolFunc = Lifetimes
        if NuseExoMolFunc > 0:
            read_hitran2exomol_path = save_path + 'conversion/HITRAN2ExoMol/'
            conversion_foldername = read_hitran2exomol_path+'/'.join(data_info)+'/'
            if os.path.exists(conversion_foldername):
                states_list = glob.glob(conversion_foldername + '__'.join(data_info[-2:]) + '.states.bz2')
                trans_list = glob.glob(conversion_foldername + '__'.join(data_info[-2:]) + '*.trans.bz2')
                if data is None:
                    hitran_df = read_hitran_parfile(
                        read_path,
                        parfile_df,
                        ConversionMinFreq,
                        ConversionMaxFreq,
                        'None',
                        'None',
                    ).reset_index(drop=True)
                else:
                    hitran_df = data.lines(
                        ConversionMinFreq,
                        ConversionMaxFreq,
                        'None',
                        'None',
                    )
                (hitran_states_col, hitran_states_fmt) = conversion_hitran2exomol(hitran_df)
            # These functions need whole states
            states_df = read_all_states(
                read_hitran2exomol_path,
                data_info,
                config.check_uncertainty,
                config.states_col,
                config.states_fmt,
            )
            if Lifetimes == 1:
                save_hitran_lifetime(read_hitran2exomol_path, states_df, hitran_states_col, hitran_states_fmt)
        # Use HITRAN functions
        # Calculate cooling functions or oscillator strengths
        NeedWholeHITRAN = PartitionFunctions + SpecificHeats + CoolingFunctions + OscillatorStrengths
        if NeedWholeHITRAN != 0:
            if data is None:
                hitran_df = read_hitran_parfile(
                    read_path,
                    parfile_df,
                    min_wn,
                    max_wn,
                    'None',
                    'None',
                ).reset_index(drop=True)
            else:
                hitran_df = data.lines(min_wn, max_wn, 'None', 'None')
        if PartitionFunctions == 1:
            save_hitran_partition_func(hitran_df, Ntemp, Tmax)
        if SpecificHeats == 1:
            save_hitran_specific_heat(hitran_df, Ntemp, Tmax)
        if CoolingFunctions == 1:
            save_hitran_cooling_func(hitran_df, Ntemp, Tmax)
        if OscillatorStrengths == 1:
            save_hitran_oscillator_strength(hitran_df)
        if (StickSpectra + CrossSections > 0):
            if data is None:
                hitran_df = read_hitran_parfile(
                    read_path,
                    parfile_df,
                    min_wn,
                    max_wn,
                    UncFilter,
                    threshold,
                ).reset_index(drop=True)
            else:
                hitran_df = data.lines(min_wn, max_wn, UncFilter, threshold)
            (hitran_linelist_df, QNs_col) = hitran_linelist_QN(hitran_df)
        if StickSpectra == 1 and CrossSections == 1:
            # Use combined function for both stick spectra and cross sections
            save_hitran_stick_spectra_cross_section(hitran_linelist_df, QNs_col, T_list, P_list, Tvib_list, Trot_list)
        elif StickSpectra == 1:
            save_hitran_stick_spectra(hitran_linelist_df, QNs_col, T_list, Tvib_list, Trot_list)
        elif CrossSections == 1:
            save_hitran_cross_section(hitran_linelist_df, T_list, P_list, Tvib_list, Trot_list)
    else:
        raise ValueError("Please add the name of the database 'ExoMol', 'ExoMolHR', 'ExoAtom', 'HITRAN', or 'HITEMP' into the input file.")
    
    print('The program total running time:')
    t_tot.end()
    print('\nFinished!')
