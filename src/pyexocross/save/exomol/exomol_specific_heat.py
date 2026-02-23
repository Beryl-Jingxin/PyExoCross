"""
Save ExoMol specific heats.

This module provides functions for calculating and saving specific heat capacities.
"""
import numpy as np
import pandas as pd
from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.base.log import print_file_info
from pyexocross.calculation.calculate_specific_heat import cal_specific_heat

# Proces Specific Heat
def save_exomol_specific_heat(states_df, Ntemp, Tmax):
    """
    Calculate and save specific heat capacities for ExoMol database.

    Computes specific heat at constant pressure at specified temperature intervals
    and saves results in .cp format file.

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

    print('Calculate specific heats.')  
    t = Timer()
    t.start()
    En = states_df['E'].astype('float').values
    gn = states_df['g'].astype('int').values
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp)) 
    specific_heat = [cal_specific_heat(En, gn, T) for T in Ts]
    t.end()
    print('Finished calculating specific heats!\n')

    print('Saving specific heats into file ...')   
    ts = Timer()    
    ts.start() 
    specific_heat_df = pd.DataFrame()
    specific_heat_df['T'] = Ts
    specific_heat_df['Specific heat'] = specific_heat 
    
    cp_folder = save_path + 'specific_heat/'
    ensure_dir(cp_folder)
    cp_path = cp_folder + '__'.join(data_info[-2:]) + '.cp'
    np.savetxt(cp_path, specific_heat_df, fmt="%8.1f %15.4f")
    ts.end()
    print_file_info('Specific heats', ['T', 'Specific heats'], ['%8.1f', '%15.4f'])
    print('Specific heats file has been saved:', cp_path, '\n')  
    print('Specific heats have been saved!\n')    
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
   