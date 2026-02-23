"""
Convert HITRAN format to ExoMol format.

This module provides functions for converting HITRAN linelist data to ExoMol
states and transitions format, including quantum number conversion.
"""
import os
import numpy as np
import pandas as pd
from ..base import Timer, ensure_dir, print_conversion_info, print_file_info
from ..calculation.calculate_para import cal_Ep
from ..process.hitran_qn import globalQNclasses, localQNgroups, separate_QN_hitran
from ..base.large_file import save_large_txt
from ..process.hitran_qn import globalQNclasses

## HITRAN to ExoMol
def convert_QNValues_hitran2exomol(hitran2exomol_states_df, GlobalQNLabel_list, LocalQNLabel_list):
    """
    Convert quantum number values from HITRAN format to ExoMol format.

    Converts symmetry labels (Gtot, Gvib, Grot) from HITRAN string format
    to ExoMol numeric codes, and converts taui from character to numeric format.

    Parameters
    ----------
    hitran2exomol_states_df : pd.DataFrame
        States DataFrame with HITRAN quantum number values
    GlobalQNLabel_list : list of str
        List of global quantum number labels
    LocalQNLabel_list : list of str
        List of local quantum number labels

    Returns
    -------
    pd.DataFrame
        States DataFrame with quantum number values converted to ExoMol format
    """
    QNLabel_list = GlobalQNLabel_list+LocalQNLabel_list
    if 'Gtot' in QNLabel_list:
        hitran2exomol_states_df["Gtot"] = (hitran2exomol_states_df["Gtot"]
                                           .str.replace('A1',"1'").str.replace('A2',"2'").str.replace("E'",'3')
                                           .str.replace('A1"','4').str.replace('A2"','5').str.replace('E"','6'))    
    if 'Gvib' in QNLabel_list:
        hitran2exomol_states_df["Gvib"] = (hitran2exomol_states_df["Gvib"]
                                           .str.replace('A1',"1'").str.replace('A2',"2'").str.replace("E'",'3')
                                           .str.replace('A1"','4').str.replace('A2"','5').str.replace('E"','6'))    
    if 'Grot' in QNLabel_list:
        hitran2exomol_states_df["Grot"] = (hitran2exomol_states_df["Grot"]
                                           .str.replace('A1',"1'").str.replace('A2',"2'").str.replace("E'",'3')
                                           .str.replace('A1"','4').str.replace('A2"','5').str.replace('E"','6'))                  
    if 'taui' in QNLabel_list:
        hitran2exomol_states_df["taui"] = hitran2exomol_states_df["taui"].str.replace('s','0').str.replace('a','1')
    return hitran2exomol_states_df

def convert_hitran2StatesTrans(hitran_df, QNu_df, QNl_df):
    """
    Convert HITRAN linelist to ExoMol states and transitions format.

    Processes upper and lower states separately, merges quantum numbers,
    deduplicates states, and creates transitions DataFrame with state IDs.

    Parameters
    ----------
    hitran_df : pd.DataFrame
        HITRAN format linelist DataFrame
    QNu_df : pd.DataFrame
        Upper state quantum numbers DataFrame
    QNl_df : pd.DataFrame
        Lower state quantum numbers DataFrame

    Returns
    -------
    tuple
        A tuple containing:
        - Jpp_df : pd.DataFrame, lower state J quantum numbers
        - hitran2exomol_states_df : pd.DataFrame, ExoMol format states DataFrame
        - hitran2exomol_trans_df : pd.DataFrame, ExoMol format transitions DataFrame
    """
    from pyexocross.core import GlobalQNLabel_list, LocalQNLabel_list
    hitran2exomol_upper_df = pd.concat([hitran_df[['A','v','gp','Unc','Ep']],QNu_df], axis=1, join='inner')
    hitran2exomol_lower_df = pd.concat([hitran_df[['A','v','gpp','Unc','Epp']],QNl_df], axis=1, join='inner')
    Jpp_df = hitran2exomol_lower_df['J']
    # hitran2exomol_upper_df['F'] = cal_F(hitran2exomol_upper_df['gp']).astype(str) 
    # hitran2exomol_lower_df['F'] = cal_F(hitran2exomol_lower_df['gpp']).astype(str) 

    hitran2exomol_upper_df.columns = list(map(lambda x: x.replace('gp','g').replace('Ep','E'), list(hitran2exomol_upper_df.columns)))
    hitran2exomol_lower_df.columns = list(map(lambda x: x.replace('gpp','g').replace('Epp','E'), list(hitran2exomol_lower_df.columns)))

    # hitran2exomol_lower_df = hitran2exomol_lower_df.reset_index(drop=True)
    # hitran2exomol_upper_df = hitran2exomol_upper_df.reset_index(drop=True)
    # if ('Sym"' in QNl_col) and ("Sym'" not in QNu_col):
    #     hitran2exomol_upper_df['Sym'] = hitran2exomol_lower_df['Sym']
    #     index_change_sym = np.where(np.array(hitran2exomol_lower_df['J']-hitran2exomol_upper_df['J'])==0.0)[0]
    #     hitran2exomol_upper_df['Sym'][index_change_sym] = (hitran2exomol_upper_df['Sym'][index_change_sym]
    #                                                        .replace('e','e2f').replace('f','e').replace('e2f','f'))
        
    hitran2exomol_ul_df = pd.concat([hitran2exomol_upper_df, hitran2exomol_lower_df], axis=0)

    hitranQNlabels = [x for x in list(hitran2exomol_ul_df.columns)[5:] if x != 'M' and x != 'm']
    hitran2exomol_states_noid = (hitran2exomol_ul_df.loc[hitran2exomol_ul_df.groupby(['g']+hitranQNlabels)['Unc'].idxmax()][['E','g','Unc']+hitranQNlabels]
                                .drop_duplicates().sort_values('E').groupby(['E','g']+hitranQNlabels)['Unc'].max().reset_index())
    hitran2exomol_states_id = hitran2exomol_states_noid
    hitran2exomol_states_id['id'] = hitran2exomol_states_noid.index+1

    states_columns_order = ['id','E','g','F','Unc']+[x for x in hitranQNlabels if x != 'F']
    hitran2exomol_states_id = hitran2exomol_states_id[states_columns_order]

    # Transitions
    upper_QNlabel = list(hitran2exomol_upper_df.columns[5:])
    upper_idAv = hitran2exomol_upper_df.merge(hitran2exomol_states_id, on=['g']+upper_QNlabel, how='inner').drop(columns=upper_QNlabel)
    upper_idAv['diffE'] = np.abs(upper_idAv['E_x']-upper_idAv['E_y'])
    upper_AvEid = upper_idAv.loc[upper_idAv.groupby(['v'])['diffE'].idxmin()][['id','A','v','E_y']].rename(columns={'id':'u','E_y':'Ep'})

    lower_QNlabel = list(hitran2exomol_lower_df.columns[5:])
    lower_idAv = hitran2exomol_lower_df.merge(hitran2exomol_states_id, on=['g']+lower_QNlabel, how='inner').drop(columns=lower_QNlabel)
    lower_idAv['diffE'] = np.abs(lower_idAv['E_x']-lower_idAv['E_y'])
    lower_AvEid = lower_idAv.loc[lower_idAv.groupby(['v'])['diffE'].idxmin()][['id','A','v','E_y',]].rename(columns={'id':'l','E_y':'Epp'})

    hitran2exomol_idAv = upper_AvEid.merge(lower_AvEid, on=['v'], how='inner').drop(columns=['A_y']).rename(columns={'A_x':'A','v':'v_x'})
    hitran2exomol_idAv['v'] = hitran2exomol_idAv['Ep']-hitran2exomol_idAv['Epp']
    diff = hitran2exomol_idAv[['u','l','A','v','v_x']].copy()
    diff.loc[:, 'diffv'] = np.abs(diff['v_x'] - diff['v'])
    hitran2exomol_trans_df = diff.loc[diff.groupby(['v'])['diffv'].idxmin()][['u','l','A','v']].sort_values('v')
    
    # States
    hitran2exomol_states_df = convert_QNValues_hitran2exomol(hitran2exomol_states_id, GlobalQNLabel_list, LocalQNLabel_list)
    # Convert HITRAN uncertainty codes to float values for ExoMol format (%12.6f)
    # HITRAN codes: 0-10 map to uncertainty values
    unc_mapping = {10: 1e-09, 9: 1e-08, 8: 1e-07, 7: 1e-06, 6: 1e-05, 
                   5: 1e-04, 4: 0.001, 3: 0.01, 2: 0.1, 1: 1.0, 0: 10.0}
    hitran2exomol_states_df['Unc'] = hitran2exomol_states_df['Unc'].map(unc_mapping).astype(float)
    
    return(Jpp_df, hitran2exomol_states_df, hitran2exomol_trans_df)

def convert_hitran2broad(hitran_df, Jpp_df):
    """
    Convert HITRAN broadening parameters to ExoMol broadening format.

    Creates air and self-broadening DataFrames with ExoMol format columns.

    Parameters
    ----------
    hitran_df : pd.DataFrame
        HITRAN format linelist DataFrame containing broadening parameters
    Jpp_df : pd.DataFrame
        Lower state J quantum numbers DataFrame

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        A tuple containing:
        - hitran2exomol_air_df : pd.DataFrame, air-broadening DataFrame
        - hitran2exomol_self_df : pd.DataFrame, self-broadening DataFrame
    """
    broad_code_df = pd.DataFrame(np.full_like(Jpp_df.astype(str),'a0'), columns=['code'])
    hitran2exomol_air_df = pd.concat([broad_code_df, hitran_df[['gamma_air','n_air']], Jpp_df], axis=1).drop_duplicates().dropna()
    hitran2exomol_self_df = pd.concat([broad_code_df, hitran_df[['gamma_self','n_air']], Jpp_df], axis=1).drop_duplicates().dropna()
    return(hitran2exomol_air_df, hitran2exomol_self_df)

def conversion_states(hitran2exomol_states_df, conversion_folder):
    """
    Save converted states data in ExoMol format.

    Formats and saves states DataFrame to .states file with appropriate
    column formats.

    Parameters
    ----------
    hitran2exomol_states_df : pd.DataFrame
        Converted states DataFrame in ExoMol format
    conversion_folder : str
        Output folder path for converted files

    Returns
    -------
    tuple of (list, list)
        A tuple containing:
        - states_col : list of str, column names
        - states_fmt : list of str, format strings for each column
    """
    from pyexocross.core import data_info, QNsformat_list, QNslabel_list
    print('Convert states data format from HITRAN to ExoMol.')  
    t = Timer()
    t.start()
    conversion_states_path = conversion_folder + '__'.join(data_info[-2:]) + '.states'
    hitranQNlabels_noF = hitran2exomol_states_df.columns[5:].tolist()
    hitranQNformats = [QNsformat_list[j] for j in [QNslabel_list.index(i) for i in hitranQNlabels_noF]]
    # states_format = ("%12s %12.6f %6s %7s %12.6f " 
    #                  + str(QNsformat_list).replace("['","").replace("']","")
    #                  .replace("'","").replace(",","").replace("d","s").replace("i","s"))
    states_format = ("%12s %12.6f %6s %7s %12.6f " 
                    + str(hitranQNformats).replace("['","").replace("']","")
                    .replace("'","").replace(",","").replace("d","s").replace("i","s"))
    np.savetxt(conversion_states_path, hitran2exomol_states_df, fmt=states_format)
    t.end()
    states_col = hitran2exomol_states_df.columns
    states_fmt = ['%12s', '%12.6f', '%6d', '%7s', '%12.6f'] + hitranQNformats
    print_file_info('Converted ExoMol states', states_col, states_fmt)
    print('Converted states file has been saved:', conversion_states_path)
    print('Converted states file has been saved!\n')   
    return(states_col, states_fmt)
    
def conversion_trans(hitran2exomol_trans_df, conversion_folder):
    """
    Save converted transitions data in ExoMol format.

    Formats and saves transitions DataFrame to compressed .trans.bz2 file.

    Parameters
    ----------
    hitran2exomol_trans_df : pd.DataFrame
        Converted transitions DataFrame in ExoMol format
    conversion_folder : str
        Output folder path for converted files
    """
    from pyexocross.core import data_info
    print('Convert transitions data format from HITRAN to ExoMol.')  
    t = Timer()
    t.start()
    conversion_trans_path = conversion_folder + '__'.join(data_info[-2:]) + '.trans.bz2'
    trans_format = "%12d %12d %10.4e %15.6f"
    save_large_txt(conversion_trans_path, hitran2exomol_trans_df, fmt=trans_format)
    t.end()
    print_file_info('Converted ExoMol transitions', ['f', 'i', 'A', 'v'], ['%12d', '%12d', '%10.4e', '%15.6f'])
    print('Converted transitions file has been saved:', conversion_trans_path)
    print('Converted transitions file has been saved!\n')  
    
def conversion_broad(hitran2exomol_air_df, hitran2exomol_self_df, conversion_folder):
    """
    Save converted broadening data in ExoMol format.

    Saves air and self-broadening DataFrames to separate .broad files
    if they contain data.

    Parameters
    ----------
    hitran2exomol_air_df : pd.DataFrame
        Air-broadening DataFrame in ExoMol format
    hitran2exomol_self_df : pd.DataFrame
        Self-broadening DataFrame in ExoMol format
    conversion_folder : str
        Output folder path for converted files
    """
    from pyexocross.core import data_info
    print('Convert broadening data format from HITRAN to ExoMol.')  
    t = Timer()
    t.start()
    nair = len(hitran2exomol_air_df)
    nself = len(hitran2exomol_self_df)
    nbroad = nair + nself
    broad_format = "%2s %6.4f %6.3f %7s"    
    if nair != 0:
        conversion_airbroad_path = conversion_folder + data_info[-2] + '__air.broad'
        np.savetxt(conversion_airbroad_path, hitran2exomol_air_df, fmt=broad_format)
        print('Converted air broadening file have been saved:', conversion_airbroad_path)
    else:
        print('No air broadening file.')
    if nself != 0:
        conversion_selfbroad_path = conversion_folder + data_info[-2] + '__self.broad'
        np.savetxt(conversion_selfbroad_path, hitran2exomol_self_df, fmt=broad_format)
        print('Converted self broadening file have been saved:', conversion_selfbroad_path)
    else:
        print('No self broadening file.')
    t.end()
    if nbroad != 0:
        hitran2exomol_broad_fmt_list = ['%2s', '%6.4f', '%6.3f', '%7s']
        print_file_info('Converted ExoMol broadening', ['Code', 'gamma_ref', 'n_L', 'J"'], hitran2exomol_broad_fmt_list)
        print('Converted broadening files have been saved!\n')  
    else:
        print('No broadening files need to be saved!\n')  

def conversion_hitran2exomol(hitran_df):
    """
    Main function to convert HITRAN database format to ExoMol format.

    Processes HITRAN linelist, converts states, transitions, and broadening
    data, and saves results in ExoMol file format (.states, .trans.bz2, .broad).

    Parameters
    ----------
    hitran_df : pd.DataFrame
        HITRAN format linelist DataFrame

    Returns
    -------
    tuple of (list, list)
        A tuple containing:
        - states_col : list of str, states column names
        - states_fmt : list of str, states format strings
    """
    from pyexocross.core import (
        data_info,
        save_path,
        ConversionMinFreq,
        ConversionMaxFreq,
        GlobalQNLabel_list,
        GlobalQNFormat_list,
        LocalQNLabel_list,
        LocalQNFormat_list,
        ConversionUnc,
        ConversionThreshold,
    )
    hitran_df = hitran_df[~hitran_df['Vp'].isin([' '*15])]
    hitran_df = hitran_df[~hitran_df['Vpp'].isin([' '*15])]
    hitran_df = hitran_df[~hitran_df['Qp'].isin([' '*15])]
    hitran_df = hitran_df[~hitran_df['Qpp'].isin([' '*15])]
    GlobalQNLabels,GlobalQNFormats = globalQNclasses(data_info)
    LocalQNupperLabels, LocalQNlowerLabels, LocalQNupperFormats, LocalQNlowerFormats = localQNgroups(data_info)
    QNu_df, QNl_df, QNu_col, QNl_col = separate_QN_hitran(hitran_df,GlobalQNLabels,LocalQNupperLabels,LocalQNlowerLabels,
                                                          GlobalQNFormats,LocalQNupperFormats,LocalQNlowerFormats)
    hitran_df['Ep'] = cal_Ep(hitran_df['Epp'].values,hitran_df['v'].values)
    
    Jpp_df, hitran2exomol_states_df, hitran2exomol_trans_df = convert_hitran2StatesTrans(hitran_df, QNu_df, QNl_df)
    
    hitran2exomol_air_df, hitran2exomol_self_df = convert_hitran2broad(hitran_df, Jpp_df)
    
    conversion_folder = save_path+'conversion/HITRAN2ExoMol/'+'/'.join(data_info)+'/' 
    ensure_dir(conversion_folder)
    print('Convert data format from HITRAN to ExoMol.')  
    print_conversion_info(ConversionMinFreq, ConversionMaxFreq, GlobalQNLabel_list, GlobalQNFormat_list, 
                          LocalQNLabel_list, LocalQNFormat_list, ConversionUnc, ConversionThreshold)
    (states_col, states_fmt) = conversion_states(hitran2exomol_states_df, conversion_folder)
    conversion_trans(hitran2exomol_trans_df, conversion_folder)
    conversion_broad(hitran2exomol_air_df, hitran2exomol_self_df, conversion_folder)
    print('Finished converting data format from HITRAN to ExoMol!\n')
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
    return(states_col, states_fmt)
    