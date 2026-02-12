import numpy as np
import pandas as pd

from src.base.utils import Timer, ensure_dir
from src.base.log import print_file_info
from src.calculation.calculate_oscillator_strength import cal_oscillator_strength
from src.plot.plot_oscillator_strength import plot_oscillator_strength


def save_hitran_oscillator_strength(hitran_df):
    """
    Main function to calculate and save oscillator strengths for HITRAN database.

    Calculates oscillator strengths from HITRAN linelist and saves results
    in .os format file.

    Parameters
    ----------
    hitran_df : pd.DataFrame
        HITRAN DataFrame with columns: A, v, gp, gpp
    """
    from pyexocross.core import save_path, data_info, database, gfORf, PlotOscillatorStrengthYN
    print('Calculate oscillator strengths.')  
    t = Timer()
    t.start()
    A = hitran_df['A'].values
    v = hitran_df['v'].values
    gp = hitran_df['gp'].values
    gpp = hitran_df['gpp'].values
    oscillator_strength_df = hitran_df[['gp', 'gpp', 'v']].copy()
    oscillator_strength_df.loc[:, 'os'] = cal_oscillator_strength(gp, gpp, A, v)
    oscillator_strength_df = oscillator_strength_df[['gp', 'gpp', 'os', 'v']]
    oscillator_strength_df = oscillator_strength_df.sort_values(by=['v'])
    t.end()
    print('Finished calculating oscillator strengths!\n')

    print('Saving oscillator strengths into file ...')      
    ts = Timer()    
    ts.start()
    os_folder = save_path + 'oscillator_strength/files/'+data_info[0]+'/'+database+'/'
    ensure_dir(os_folder)
    os_path = os_folder+'__'.join(data_info[-2:])+'_'+gfORf.lower()+'.os' 
    os_format = "%7.1f %7.1f %10.4E %15.6f"
    np.savetxt(os_path, oscillator_strength_df, fmt=os_format)
    ts.end()
    print_file_info('Oscillator strengths', oscillator_strength_df.columns, ['%7.1f', '%7.1f', '%10.4E', '%15.6f'])
    print('Wavenumber in unit of cm⁻¹')
    print('Oscillator strengths file has been saved:', os_path, '\n')  
    
    # Plot oscillator strengths and save it as .png.
    if PlotOscillatorStrengthYN == 'Y':
        plot_oscillator_strength(oscillator_strength_df)  
    print('\nOscillator strengths have been saved!\n') 
    print('* * * * * - - - - - * * * * * - - - - - * * * * * - - - - - * * * * *\n')
