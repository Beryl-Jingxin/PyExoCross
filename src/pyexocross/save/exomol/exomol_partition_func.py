"""
Save ExoMol partition functions.

This module provides functions for calculating and saving partition functions.
"""
import numpy as np
import pandas as pd
from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.base.log import print_file_info
from pyexocross.calculation.calculate_partition_func import cal_partition_func

# Process Partition Function
def save_exomol_partition_func(states_df, Ntemp, Tmax):
    """
    Calculate and save partition functions for ExoMol database.

    Computes LTE partition functions at specified temperature intervals
    and saves results in .pf format file.

    Parameters
    ----------
    states_df : pd.DataFrame
        Complete states DataFrame with 'E' and 'g' columns
    Ntemp : int
        Temperature step interval
    Tmax : int
        Maximum temperature in Kelvin
    """
    # Import legacy-style configuration variables from core (set via Config.to_globals()).
    from pyexocross.core import save_path, data_info  

    print('Calculate partition functions.')  
    t = Timer()
    t.start()
    En = states_df['E'].astype('float').values
    gn = states_df['g'].astype('int').values
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp)) 
    partition_func = [cal_partition_func(En, gn, T) for T in Ts]
    t.end()
    print('Finished calculating partition functions!\n')

    print('Saving partition functions into file ...')   
    ts = Timer()    
    ts.start()     
    partition_func_df = pd.DataFrame()
    partition_func_df['T'] = Ts
    partition_func_df['Partition function'] = partition_func
    
    pf_folder = save_path + 'partition/'
    ensure_dir(pf_folder)
    pf_path = pf_folder + '__'.join(data_info[-2:]) + '.pf'
    np.savetxt(pf_path, partition_func_df, fmt="%8.1f %15.4f")
    ts.end()
    print_file_info('Partition functions', ['T', 'Partition function'], ['%8.1f', '%15.4f'])
    print('Partition functions file has been saved:', pf_path, '\n') 
    print('Partition functions have been saved!\n')  
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')

