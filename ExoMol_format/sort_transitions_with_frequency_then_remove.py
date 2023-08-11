# Import all what we need.
# encoding: utf-8
import os
import bz2
import csv
import glob
import requests
import numpy as np
import pandas as pd
import numexpr as ne
from io import StringIO
from tqdm.notebook import tqdm
import warnings
warnings.simplefilter("ignore", np.ComplexWarning)
pd.options.mode.chained_assignment = None


molecule = 'ZrO'
isotopologue = '96Zr-16O'
dataset = 'ZorrO'
ncolumn = 3
read_path = '/mnt/data/exomol/exomol3_data/'


def read_all_states(read_path):
    print('Reading states ...')
    s_df = dict()
    states_df = pd.DataFrame()
    states_filenames = glob.glob(read_path + molecule + '/' + isotopologue + '/' + dataset 
                                 + '/' + isotopologue + '__' + dataset + '.states.bz2')
    for states_filename in states_filenames:
        s_df[states_filename] = pd.read_csv(states_filename, compression='bz2', sep='\s+', header=None,
                                            chunksize=1_000_000, iterator=True, low_memory=False, dtype=object)
        for chunk in s_df[states_filename]:
            states_df = pd.concat([states_df, chunk])   
    states_df = states_df.rename(columns={0:'id',1:'E'})  
    
    print('Finished reading states!\n')                       
    return(states_df)

def get_transfiles(read_path):
    # Get all the transitions files from the folder including the older version files which are named by vn(version number).
    trans_filepaths_all = glob.glob(read_path + molecule + '/' + isotopologue + '/' + dataset + '/' + '*trans.bz2')
    num_transfiles_all = len(trans_filepaths_all)    # The number of all transitions files including the older version files.
    trans_filepaths = []                             # The list of the lastest transitions files.
    for i in range(num_transfiles_all):
        split_version = trans_filepaths_all[i].split('__')[-1].split('.')[0].split('_')    # Split the filenames.
        num = len(split_version)
        if num == 1:     
            if split_version[0] == dataset:        
                trans_filepaths.append(trans_filepaths_all[i])
            if len(split_version[0].split('-')) == 2:
                trans_filepaths.append(trans_filepaths_all[i])
    return(trans_filepaths)    

def read_all_trans(read_path):
    print('Reading all transitions ...')
    t_df = dict()
    all_trans_df = pd.DataFrame()
    trans_filepaths = get_transfiles(read_path)
    for trans_filename in tqdm(trans_filepaths, position=0, leave=True, desc='Reading transitions'):
        t_df[trans_filename] = pd.read_csv(trans_filename, compression='bz2', sep='\s+', header=None,
                                           chunksize=1_000_000, iterator=True, low_memory=False)
        for chunk in t_df[trans_filename]:
            all_trans_df = pd.concat([all_trans_df,chunk])
    ncolumn = len(all_trans_df.columns)
    if ncolumn == 3: 
        trans_col_name={0:'u', 1:'l', 2:'A'}
    else:
        trans_col_name={0:'u', 1:'l', 2:'A', 3:'v'}
    all_trans_df = all_trans_df.rename(columns=trans_col_name)                
    print('Finished reading all transitions!\n')                         
    return(all_trans_df, ncolumn)

def cal_v(Ep, Epp):
    v = ne.evaluate('abs(Ep - Epp)')
    return(v)

states_df = read_all_states(read_path)
(all_trans_df, ncolumn) = read_all_trans(read_path)
    
id_u = all_trans_df['u'].values
id_l = all_trans_df['l'].values
states_df['id'] = pd.to_numeric(states_df['id'])
states_df.set_index(['id'], inplace=True, drop=False)
id_s = states_df['id']
all_trans_df.set_index(['u'], inplace=True, drop=False)
id_us = list(set(id_u).intersection(set(id_s)))
trans_us_df = all_trans_df.loc[id_us]
id_l = trans_us_df['l'].values
id_ls = list(set(id_l).intersection(set(id_s)))
trans_us_df.set_index(['l'], inplace=True, drop=False)
trans_s_df = trans_us_df.loc[id_ls]
id_su = trans_s_df['u'].values
id_sl = trans_s_df['l'].values
states_u_df = states_df.loc[id_su]
states_l_df = states_df.loc[id_sl]

Ep = states_u_df['E'].values.astype('float')
A = trans_s_df['A'].values.astype('float')

if ncolumn == 4:
    v = trans_s_df['v'].values.astype('float')
    trans_format = "%12d %12d %10.4e %15.6f"
else:
    Epp = states_l_df['E'].astype('float')     # Upper state energy
    v = cal_v(Ep, Epp) 
    trans_format = "%12d %12d %10.4e"

trans_s_df['v'] = v
trans_s_df.sort_values('v',inplace=True)
trans_s_df3 = trans_s_df[['u','l','A']]

trans_path = (read_path + molecule + '/' + isotopologue + '/' + dataset + '/' + 
              isotopologue + '__' + dataset + '.trans.bz2')
np.savetxt(trans_path, trans_s_df3, fmt=trans_format)
