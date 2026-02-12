import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

from src.base.utils import Timer, ensure_dir
from src.base.log import log_tqdm, print_file_info
from src.calculation.calculate_para import cal_Ep
from src.calculation.calculate_cooling_func import cal_cooling_func
from src.database.load_hitran import read_hitran_pf


def save_hitran_cooling_func(hitran_df, Ntemp, Tmax):
    """
    Main function to calculate and save cooling functions for HITRAN database.

    Calculates cooling functions at specified temperature intervals from HITRAN
    linelist and saves results in .cf format file.

    Parameters
    ----------
    hitran_df : pd.DataFrame
        HITRAN DataFrame with columns: A, v, gp, Epp
    Ntemp : int
        Temperature step interval
    Tmax : int
        Maximum temperature in Kelvin
    """
    from pyexocross.core import ncputrans, save_path, data_info
    print('Calculate cooling functions.')  
    t = Timer()
    t.start()
    A = hitran_df['A'].values
    v = hitran_df['v'].values
    gp = hitran_df['gp'].values
    Ep = cal_Ep(hitran_df['Epp'].values,v)
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp)) 
    Qs = read_hitran_pf(Ts)
    # Parallelize per-temperature cooling evaluation
    with ProcessPoolExecutor(max_workers=ncputrans) as executor:
        futures = [executor.submit(cal_cooling_func, A, v, Ep, gp, T_val, Q_val)
                   for T_val, Q_val in log_tqdm(zip(Ts, Qs), desc='Calculating')]
        cooling_func = [future.result() for future in futures]
    t.end()
    print('Finished calculating cooling functions!\n')

    print('Saving cooling functions into file ...')   
    ts = Timer()    
    ts.start()      
    cooling_func_df = pd.DataFrame()
    cooling_func_df['T'] = Ts
    cooling_func_df['cooling function'] = cooling_func

    cf_folder = save_path + 'cooling/'
    ensure_dir(cf_folder)
    cf_path = cf_folder + '__'.join(data_info[-2:]) + '.cf' 
    np.savetxt(cf_path, cooling_func_df, fmt="%8.1f %20.8E")
    ts.end()
    print_file_info('Cooling functions', ['T', 'Cooling function'], ['%8.1f','%20.8E'])
    print('Cooling functions file has been saved:', cf_path, '\n')
    print('Cooling functions have been saved!\n')     
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
