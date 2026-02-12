"""
Process HITRAN quantum numbers.

This module provides functions for processing and converting HITRAN quantum
number formats.
"""
import numpy as np
import pandas as pd
import src.base.input as _input
from src.calculation.calculate_para import cal_Ep_hitran
from src.process.filter_qn import QNfilter_linelist

## Process HITRAN quantum numbers
# Global quantum numbers
def globalQNclasses(data_info):
    """
    Determine global quantum number labels and formats for HITRAN molecule.

    Maps molecule/isotopologue to appropriate global quantum number class
    based on HITRAN classification system.

    Parameters
    ----------
    data_info : list or tuple
        Data information containing [molecule, isotopologue]

    Returns
    -------
    tuple of (list, list)
        A tuple containing:
        - GlobalQNLabels : list of str, global quantum number labels
        - GlobalQNFormats : list of str, format strings for each label
    """
    globalQNclass1a = {'class':['CO','HF','HCl','HBr','HI','N2','NO+','NO_p','H2','CS'],
                       'label': ['none','v1'],
                       'format':['%13s','%2d']}
    globalQNclass1b = {'class':['O2','NO','OH','ClO','SO','S2'],
                       'label':['none','X','Omega','none','v1'],
                       'format':['%6s','%2s','%3s','%2s','%2d']}
    globalQNclass2a = {'class':['CO2'],
                       'label':['none','v1','v2','l2','v3','r'],
                       'format':['%6s','%2d','%2d','%2d','%2d','%1d']}
    globalQNclass2b = {'class':['N2O','OCS','HCN','CS2'],
                       'label':['none','v1','v2','l2','v3'],
                       'format':['%7s','%2d','%2d','%2d','%2d']}
    globalQNclass3  = {'class':['H2O','O3','SO2','NO2','HOCl','H2S','HO2','HOBr'],
                       'label':['none','v1','v2','v3'],
                       'format':['%9s','%2d','%2d','%2d']}
    globalQNclass4a = {'class':['PH3','NF3'],
                       'label':['none','v1','v2','v3','v4','S'],
                       'format':['%5s','%2d','%2d','%2d','%2d','%2s']}
    globalQNclass4b = {'class':['14N-1H3','15N-1H3'],
                       'label':['none','v1','v2','v3','v4','none','l3','l4','none','l','none','Gvib'],
                       'format':['%1s','%1d','%1d','%1d','%1d','%1s','%1d','%1d','%1s','%1d','%1s','%4s']}
    globalQNclass5a = {'class':['C2H2'],
                       'label':['none','v1','v2','v3','v4','v5','l4','l5','+/-','none','S'],
                       'format':['%1s','%1d','%1d','%1d','%2d','%2d','%2d','%2d','%1s','%1s','%1s']}
    globalQNclass5b = {'class':['C4H2'],
                       'label':['none','v1','v2','v3','v4','v5','v6','v7','v8','v9','none','Sym','none','S'],
                       'format':['%1s','%1d','%1d','%1d','%1d','%1d','%1d','%1d','%1d','%1d','%1s','%1s','%1s','%2s']}
    globalQNclass5c = {'class':['HC3N'],
                       'label':['none','v1','v2','v3','v4','v5','v6','v7','l5','l6','l7'],
                       'format':['%2s','%1d','%1d','%1d','%1d','%1d','%1d','%1d','%2d','%2d','%2d']}
    globalQNclass5d = {'class':['C2N2'],
                       'label':['v1','v2','v3','v4','v5','l','+/-','r','S'],
                       'format':['%2d','%2d','%2d','%2d','%2d','%2d','%1s','%1d','%1s']}
    globalQNclass6a = {'class':['H2CO','COF2','COCl2'],
                       'label':['none','v1','v2','v3','v4','v5','v6'],
                       'format':['%3s','%2d','%2d','%2d','%2d','%2d','%2d']}
    globalQNclass6b = {'class':['H2O2'],
                       'label':['none','v1','v2','v3','n','Tau','v5','v6'],
                       'format':['%3s','%2d','%2d','%2d','%1d','%1d','%2d','%2d']}
    globalQNclass7a = {'class':['SO3','CH3'],
                       'label':['v1','v2','v3','l3','v4','l4','Gvib'],
                       'format':['%2d','%2d','%2d','%2d','%2d','%2d','%3s']}
    globalQNclass7b = {'class':['H3+','H3_p'],
                       'label':['none','v1','v2','l2','none','Gtot'],
                       'format':['%5s','%2d','%2d','%2d','%1s','%3s']}
    globalQNclass8  = {'class':['12C-1H4','13C-1H4','CF4','GeH4'],
                       'label':['none','v1','v2','v3','v4','n','C'],
                       'format':['%3s','%2d','%2d','%2d','%2d','%2s','%2s']}
    globalQNclass9  = {'class':['12C-1H3-1D','12C-1H3-2H','13C-1H3-1D','13C-1H3-2H','HNO3','CH3Cl','C2H6','SF6','HCOOH','ClONO2','C2H4','CH3OH','CH3Br','CH3CN','CH3F','CH3I','COFCl','HONO','ClNO2'],
                       'label':['VibBand'],
                       'format':['%15s']}
    globalQNclass = [globalQNclass1a,globalQNclass1b,globalQNclass2a,globalQNclass2b,globalQNclass3,
                     globalQNclass4a,globalQNclass4b,globalQNclass5a,globalQNclass5b,globalQNclass5c,globalQNclass5d,
                     globalQNclass6a,globalQNclass6b,globalQNclass7a,globalQNclass7b,globalQNclass8,globalQNclass9]
    count = 0
    n_classes = len(globalQNclass) 
    molecule = data_info[0]
    isotopologue = data_info[1]
    for gQNclass in globalQNclass:
        if molecule in gQNclass.get('class'):
            GlobalQNLabels = gQNclass.get('label')
            GlobalQNFormats = gQNclass.get('format')   
        elif isotopologue in gQNclass.get('class'):
            GlobalQNLabels = gQNclass.get('label')
            GlobalQNFormats = gQNclass.get('format')
        else:
            count += 1
    if count ==  n_classes:
        GlobalQNLabels = _input.GlobalQNLabel_list
        GlobalQNFormats = _input.GlobalQNFormat_list
    return(GlobalQNLabels,GlobalQNFormats)

# Local quantum numbers
def localQNgroups(data_info):
    """
    Determine local quantum number labels and formats for HITRAN molecule.

    Maps molecule/isotopologue to appropriate local quantum number group
    for upper and lower states based on HITRAN classification system.

    Parameters
    ----------
    data_info : list or tuple
        Data information containing [molecule, isotopologue]

    Returns
    -------
    tuple of (list, list, list, list)
        A tuple containing:
        - LocalQNupperLabels : list of str, upper state local quantum number labels
        - LocalQNlowerLabels : list of str, lower state local quantum number labels
        - LocalQNupperFormats : list of str, format strings for upper state labels
        - LocalQNlowerFormats : list of str, format strings for lower state labels
    """
    localQNgroup1  = {'group':['H2O','O3','SO2','NO2','HNO3','H2CO','HOCl','H2O2','COF2','H2S','HCOOH','HO2','ClONO2','HOBr','C2H4','COCl2','COFCl','HONO','ClNO2'],
                      'ulabel': ['J','Ka','Kc','F','Sym'],
                      'uformat':['%3d','%3d','%3d','%5s','%1s'],
                      'llabel': ['J','Ka','Kc','F','Sym'],
                      'lformat':['%3d','%3d','%3d','%5s','%1s']}
    localQNgroup2a = {'group':['CO2','N2O','CO','HF','HCl','HBr','HI','OCS','N2','HCN','C2H2','NO+','NO_p','HC3N','H2','CS','C2N2','CS2'],
                      'ulabel':['m','none','F'],
                      'uformat':['%1s','%9s','%5s'],
                      'llabel': ['none','Br','J','Sym','F'],
                      'lformat':['%5s','%1s','%3d','%1s','%5s']}
    localQNgroup2b = {'group':['C4H2'],
                      'ulabel':['l6','l7','l8','l9','none'],
                      'uformat':['%2s','%2s','%2s','%2s','%7s'],
                      'llabel': ['l6','l7','l8','l9','none','Br','J','Sym','none'],
                      'lformat':['%2s','%2s','%2s','%2s','%1s','%1s','%3d','%1s','%1s']}
    localQNgroup3  = {'group':['12C-1H4','13C-1H4','SF6','CF4','GeH4'],
                      'ulabel':['none','J','C','alpha','F'],
                      'uformat':['%2s','%3d','%2s','%3d','%5s'],
                      'llabel': ['none','J','C','alpha','F'],
                      'lformat':['%2s','%3d','%2s','%3d','%5s']}
    localQNgroup4a = {'group':['12C-1H3-1D','12C-1H3-2H','13C-1H3-1D','13C-1H3-2H','CH3Cl','PH3','CH3OH','CH3Br','CH3CN','CH3F','CH3I','NF3'],
                      'ulabel':['J','K','l','C','Sym','F'],
                      'uformat':['%3d','%3d','%2d','%2s','%1s','%4s'],
                      'llabel': ['J','K','l','C','Sym','F'],
                      'lformat':['%3d','%3d','%2d','%2s','%1s','%4s']}
    localQNgroup4b = {'group':['14N-1H3','15N-1H3'],
                      'ulabel':['J','K','l','none','Grot','Gtot','none'],
                      'uformat':['%2d','%3d','%2d','%1s','%3s','%3s','%1s'],
                      'llabel': ['J','K','l','none','Grot','Gtot','none'],
                      'lformat':['%2d','%3d','%2d','%1s','%3s','%3s','%1s']}
    localQNgroup4c = {'group':['C2H6'],
                      'ulabel':['J','K','l','Sym','F'],
                      'uformat':['%3d','%3d','%2d','%3s','%4s'],
                      'llabel': ['J','K','l','Sym','F'],
                      'lformat':['%3d','%3d','%2d','%3s','%4s']}
    localQNgroup5a = {'group':['SO3','CH3'],
                      'ulabel':['J','K','none','Grot','Gtot','none'],
                      'uformat':['%3d','%3d','%2s','%3s','%3s','%1s'],
                      'llabel': ['J','K','none','Grot','Gtot','none'],
                      'lformat':['%3d','%3d','%2s','%3s','%3s','%1s']}
    localQNgroup5b = {'group':['H3+','H3_p'],
                      'ulabel':['J','Sym','G','U','K','none'],
                      'uformat':['%3d','%2s','%3d','%2s','%3s','%2s'],
                      'llabel': ['J','Sym','G','U','K','none'],
                      'lformat':['%3d','%2s','%3d','%2s','%3s','%2s']}
    localQNgroup6  = {'group':['O2','SO','S2'],
                      'ulabel':['none','F'],
                      'uformat':['%10s','%5s'],
                      'llabel': ['none','Br1','N','Br2','J','F','M'],
                      'lformat':['%1s','%1s','%3d','%1s','%3d','%5s','%1s']}
    localQNgroup7a = {'group':['NO','ClO'],
                      'ulabel':['none','none','F'],     # m, none, F
                      'uformat':['%1s','%9s','%5s'],
                      'llabel': ['none','Br','J','Sym','F'],
                      'lformat':['%2s','%2s','%5.1f','%1s','%5s']}
    localQNgroup7b = {'group':['OH'],
                      'ulabel':['none','F'],
                      'uformat':['%10s','%5s'],
                      'llabel': ['none','Br','J','Sym','F'],
                      'lformat':['%1s','%2s','%5.1f','%2s','%5s']}
    localQNgroup = [localQNgroup1,localQNgroup2a,localQNgroup2b,localQNgroup3,
                    localQNgroup4a,localQNgroup4b,localQNgroup4c,localQNgroup5a,
                    localQNgroup5b,localQNgroup6,localQNgroup7a,localQNgroup7b]
    count = 0
    n_groups = len(localQNgroup)
    molecule = data_info[0]
    isotopologue = data_info[1]
    for lQNgroup in localQNgroup:
        if molecule in lQNgroup.get('group'):
            LocalQNupperLabels = lQNgroup.get('ulabel')
            LocalQNupperFormats = lQNgroup.get('uformat')
            LocalQNlowerLabels = lQNgroup.get('llabel')
            LocalQNlowerFormats = lQNgroup.get('lformat')
        elif isotopologue in lQNgroup.get('group'):
            LocalQNupperLabels = lQNgroup.get('ulabel')
            LocalQNupperFormats = lQNgroup.get('uformat')
            LocalQNlowerLabels = lQNgroup.get('llabel')
            LocalQNlowerFormats = lQNgroup.get('lformat')
        else:
            count += 1
    if count == n_groups:
        LocalQNupperLabels = _input.LocalQNLabel_list
        LocalQNupperFormats = _input.LocalQNFormat_list
        LocalQNlowerLabels = _input.LocalQNLabel_list
        LocalQNlowerFormats = _input.LocalQNFormat_list     
    return(LocalQNupperLabels, LocalQNlowerLabels, LocalQNupperFormats, LocalQNlowerFormats)

def separate_QN_hitran(hitran_df,GlobalQNLabels,LocalQNupperLabels,LocalQNlowerLabels,
                       GlobalQNFormats,LocalQNupperFormats,LocalQNlowerFormats):
    """
    Separate global and local quantum numbers from HITRAN fixed-format strings.

    Parses V'V" and Q'Q" columns from HITRAN DataFrame according to format
    specifications and extracts individual quantum number values.

    Parameters
    ----------
    hitran_df : pd.DataFrame
        HITRAN DataFrame with Vp, Vpp, Qp, Qpp columns
    GlobalQNLabels : list of str
        List of global quantum number labels
    LocalQNupperLabels : list of str
        List of upper state local quantum number labels
    LocalQNlowerLabels : list of str
        List of lower state local quantum number labels
    GlobalQNFormats : list of str
        Format strings for global quantum numbers
    LocalQNupperFormats : list of str
        Format strings for upper state local quantum numbers
    LocalQNlowerFormats : list of str
        Format strings for lower state local quantum numbers

    Returns
    -------
    tuple
        A tuple containing:
        - QNu_df : pd.DataFrame, upper state quantum numbers
        - QNl_df : pd.DataFrame, lower state quantum numbers
        - QNu_col : list of str, upper state column names with ' suffix
        - QNl_col : list of str, lower state column names with " suffix
    """
    GQN_format = [int(s) for s in str(GlobalQNFormats).replace('.1','') if s.isdigit()]
    LQNu_format = [int(s) for s in str(LocalQNupperFormats).replace('.1','') if s.isdigit()]
    LQNl_format = [int(s) for s in str(LocalQNlowerFormats).replace('.1','') if s.isdigit()]
    hitran_df['Qp'] = [x.replace(' .5','0.5') for x in hitran_df['Qp'].values]
    hitran_df['Qpp'] = [x.replace(' .5','0.5') for x in hitran_df['Qpp'].values]
    hitran_df['Vp'] = [x.replace(' .5','0.5') for x in hitran_df['Vp'].values]
    hitran_df['Vpp'] = [x.replace(' .5','0.5') for x in hitran_df['Vpp'].values]
    # Global quantum numbers    
    n_GQN = len(GlobalQNLabels)
    GQNlsum = []
    GQNrsum = []
    for i in range(n_GQN):
        GQNlsum.append(sum(GQN_format[:i]))
        GQNrsum.append(sum(GQN_format[:i+1]))  
    GQNu_df = pd.DataFrame(columns=GlobalQNLabels)  
    GQNl_df = pd.DataFrame(columns=GlobalQNLabels)  
    for i in range(n_GQN):
        GQNu_df[GlobalQNLabels[i]] = hitran_df['Vp'].map(lambda x: x[GQNlsum[i]:GQNrsum[i]]) 
        GQNl_df[GlobalQNLabels[i]] = hitran_df['Vpp'].map(lambda x: x[GQNlsum[i]:GQNrsum[i]]) 
    if 'none' in GlobalQNLabels:
        GQNu_df = GQNu_df.drop(columns=['none'])
        GQNl_df = GQNl_df.drop(columns=['none'])

    # Local quantum numbers    
    n_LQNu = len(LocalQNupperLabels)
    LQNulsum = []
    LQNursum = []
    for i in range(n_LQNu):
        LQNulsum.append(sum(LQNu_format[:i]))
        LQNursum.append(sum(LQNu_format[:i+1]))
    n_LQNl = len(LocalQNlowerLabels)
    LQNllsum = []
    LQNlrsum = []
    for i in range(n_LQNl):
        LQNllsum.append(sum(LQNl_format[:i]))
        LQNlrsum.append(sum(LQNl_format[:i+1]))
    LQNu_df = pd.DataFrame(columns=LocalQNupperLabels)  
    for i in range(n_LQNu):
        LQNu_df[LocalQNupperLabels[i]] = hitran_df['Qp'].map(lambda x: x[LQNulsum[i]:LQNursum[i]]) 
    if 'none' in LocalQNupperLabels:
        LQNu_df = LQNu_df.drop(columns=['none'])
    LQNl_df = pd.DataFrame(columns=LocalQNlowerLabels)  
    for i in range(n_LQNl):
        LQNl_df[LocalQNlowerLabels[i]] = hitran_df['Qpp'].map(lambda x: x[LQNllsum[i]:LQNlrsum[i]]) 
    if 'none' in LocalQNlowerLabels:
        LQNl_df = LQNl_df.drop(columns=['none'])
        
    if 'Br' in LocalQNupperLabels:
        LQNu_df = LQNu_df.drop(columns=['Br'])    
    if 'Br' in LocalQNlowerLabels:
        # LQNl_df['Br'] = LQNl_df['Br'].map(lambda x: x[1]).replace(['Q','P','R','O','S'],[0,-1,1,-2,2])
        br_char = LQNl_df['Br'].map(lambda x: x[1])
        br_map = {'Q': 0, 'P': -1, 'R': 1, 'O': -2, 'S': 2}
        LQNl_df['Br'] = pd.to_numeric(br_char.map(br_map), errors='coerce')
        LQNl_df['J'] = pd.to_numeric(LQNl_df['J'])
        LQNu_df['J'] = LQNl_df['Br'] + LQNl_df['J']
        LQNl_df = LQNl_df.drop(columns=['Br'])

    if 'F' not in LocalQNupperLabels:
        LQNu_df['F'] = LQNu_df['J']
        LQNl_df['F'] = LQNl_df['J']
        
    QNu_df = pd.concat([GQNu_df, LQNu_df], axis='columns')
    QNl_df = pd.concat([GQNl_df, LQNl_df], axis='columns')
    QNu_col = [i+"'" for i in list(QNu_df.columns)] 
    QNl_col = [i+'"' for i in list(QNl_df.columns)] 
    return(QNu_df, QNl_df, QNu_col, QNl_col)

def hitran_linelist_QN(hitran_df):
    """
    Process HITRAN linelist to extract and format quantum numbers.

    Determines quantum number classes, separates quantum numbers from
    fixed-format strings, calculates upper state energies, and optionally
    applies quantum number filtering.

    Parameters
    ----------
    hitran_df : pd.DataFrame
        HITRAN DataFrame with Vp, Vpp, Qp, Qpp columns

    Returns
    -------
    tuple of (pd.DataFrame, list)
        A tuple containing:
        - hitran_linelist_df : pd.DataFrame, processed linelist with quantum numbers
        - QNs_col : list of str, quantum number column names used
    """
    from pyexocross.core import data_info, QNsFilter, QNs_label, QNs_value
    GlobalQNLabels,GlobalQNFormats = globalQNclasses(data_info)
    LocalQNupperLabels, LocalQNlowerLabels, LocalQNupperFormats, LocalQNlowerFormats = localQNgroups(data_info)
    QNu_df, QNl_df, QNu_col, QNl_col = separate_QN_hitran(hitran_df,GlobalQNLabels,LocalQNupperLabels,LocalQNlowerLabels,
                                                          GlobalQNFormats,LocalQNupperFormats,LocalQNlowerFormats)
    Ep_df = cal_Ep_hitran(hitran_df)
    hitran_linelist_df = hitran_df[['v','S','A','gamma_air','gamma_self','Epp','n_air','delta_air','gp','gpp']]
    hitran_linelist_df = pd.concat([hitran_linelist_df, Ep_df, QNu_df, QNl_df], axis='columns')
    hitran_main_colname = ['v','S','A','gamma_air','gamma_self','E"','n_air','delta_air',"g'",'g"',"E'"]
    QN_col = QNu_col + QNl_col
    hitran_linelist_colname = hitran_main_colname + QN_col
    hitran_linelist_df.columns = hitran_linelist_colname
    # Do quantum number filter.
    if QNsFilter !=[]:   
        if 'J' not in QNs_label:
            QNs_col = [i+"'" for i in QNs_label] + ["J'"] + [i+'"' for i in QNs_label] + ['J"']
        else:
            QNs_col = [i+"'" for i in QNs_label] + [i+'"' for i in QNs_label]
            hitran_linelist_df = hitran_linelist_df[hitran_main_colname + QNs_col]
        hitran_linelist_df = QNfilter_linelist(hitran_linelist_df, QNs_value, QNs_label)
    else:
        QNs_col = []  
    return(hitran_linelist_df, QNs_col) 
    