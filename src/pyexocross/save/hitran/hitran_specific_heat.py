import numpy as np
import pandas as pd

from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.base.log import print_file_info
from pyexocross.calculation.calculate_specific_heat import cal_specific_heat
from pyexocross.process.hitran_states import hitran_state_arrays


def save_hitran_specific_heat(hitran_df, Ntemp, Tmax):
    """
    Calculate and save specific heat capacities for HITRAN database.

    HITRAN .par files are transition records rather than complete state lists,
    so this uses the same approximate state reconstruction as HITRAN partition
    function calculations.
    """
    from pyexocross.core import save_path, data_info

    print('Calculate specific heats.')
    t = Timer()
    t.start()
    En, gn = hitran_state_arrays(hitran_df)
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp))
    specific_heat = np.array([cal_specific_heat(En, gn, T) for T in Ts], dtype=float)
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
