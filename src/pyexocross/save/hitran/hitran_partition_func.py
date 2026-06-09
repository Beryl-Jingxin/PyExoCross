import numpy as np
import pandas as pd
from pyexocross.base.utils import Timer, ensure_dir
from pyexocross.base.log import print_file_info
from pyexocross.calculation.calculate_partition_func import cal_partition_func


def normalize_state_symmetry(value):
    """
    Normalize HITRAN lower-state symmetry labels for state grouping.

    Some HITRAN class-1 local QN fields store two parity/symmetry characters
    such as "ef" or "fe". For a lower-state partition-function key, the second
    character is the lower-state symmetry; single-character labels are already
    state labels and are left unchanged.
    """
    text = str(value).strip()
    if len(text) == 2 and all(char in 'ef+-' for char in text):
        return text[1]
    return text


def build_hitran_partition_states(hitran_linelist_df):
    """
    Build approximate HITRAN states for partition-function calculations.

    HITRAN .par files are transition records, not state lists. Hyperfine
    components and multiple branches can make one physical lower state appear
    many times. Collapse records by lower-state quantum numbers excluding F,
    use the lowest listed component energy, and sum the unique F statistical
    weights for that physical state.
    """
    f_col = 'F"' if 'F"' in hitran_linelist_df.columns else None
    qn_cols = [
        col for col in hitran_linelist_df.columns
        if col.endswith('"') and col not in ('E"', 'g"', 'F"')
    ]
    if not qn_cols:
        states_df = hitran_linelist_df[['E"', 'g"']].rename(columns={'E"': 'E', 'g"': 'g'}).copy()
        states_df['E'] = states_df['E'].astype('float')
        states_df['g'] = states_df['g'].astype('int')
        return states_df.groupby('E', as_index=False)['g'].max()

    value_cols = ['E"', 'g"'] + qn_cols
    if f_col is not None:
        value_cols.append(f_col)
    records = (
        hitran_linelist_df[value_cols]
        .rename(columns={'E"': 'E', 'g"': 'g'})
        .copy()
    )
    records['E'] = records['E'].astype('float')
    records['g'] = records['g'].astype('int')
    if 'Sym"' in records.columns:
        records['Sym"'] = records['Sym"'].map(normalize_state_symmetry)

    rows = []
    for _, group in records.groupby(qn_cols, dropna=False):
        if f_col is not None:
            blank_f = group[f_col].astype(str).str.strip() == ''
            if blank_f.any():
                g = group.loc[blank_f, 'g'].max()
            else:
                g = group['g'].drop_duplicates().sum()
        else:
            g = group['g'].max()
        rows.append({
            'E': group['E'].min(),
            'g': g,
        })
    states_df = pd.DataFrame(rows)

    # Very high states can appear only as upper states in finite HITRAN .par
    # files. Add only upper groups whose comparable lower-state QNs are absent;
    # this avoids double-counting normal lower-state coverage.
    comparable_lower_cols = [
        col for col in qn_cols
        if col != 'Sym"' and col[:-1] + "'" in hitran_linelist_df.columns
    ]
    if comparable_lower_cols:
        lower_keys = set(
            tuple(str(row[col]).strip() for col in comparable_lower_cols)
            for _, row in hitran_linelist_df[comparable_lower_cols].drop_duplicates().iterrows()
        )
        upper_group_cols = [col[:-1] + "'" for col in comparable_lower_cols]
        upper_f_col = "F'" if "F'" in hitran_linelist_df.columns else None
        upper_value_cols = ["E'", "g'"] + upper_group_cols
        if upper_f_col is not None:
            upper_value_cols.append(upper_f_col)
        upper_records = (
            hitran_linelist_df[upper_value_cols]
            .rename(columns={"E'": 'E', "g'": 'g'})
            .copy()
        )
        upper_records['E'] = upper_records['E'].astype('float')
        upper_records['g'] = upper_records['g'].astype('int')

        upper_rows = []
        for key, group in upper_records.groupby(upper_group_cols, dropna=False):
            if not isinstance(key, tuple):
                key = (key,)
            normalized_key = tuple(str(value).strip() for value in key)
            if normalized_key in lower_keys:
                continue
            if upper_f_col is not None:
                blank_f = group[upper_f_col].astype(str).str.strip() == ''
                if blank_f.any():
                    g = group.loc[blank_f, 'g'].max()
                else:
                    g = group['g'].drop_duplicates().sum()
            else:
                g = group['g'].max()
            upper_rows.append({
                'E': group['E'].min(),
                'g': g,
            })
        if upper_rows:
            states_df = pd.concat([states_df, pd.DataFrame(upper_rows)], ignore_index=True)

    return states_df


def calculate_linelist_partition_func(hitran_df, Ts):
    from pyexocross.process.hitran_qn import hitran_linelist_QN

    hitran_linelist_df, _ = hitran_linelist_QN(hitran_df)
    states_df = build_hitran_partition_states(hitran_linelist_df)

    En = states_df['E'].astype('float').values
    gn = states_df['g'].astype('int').values
    return np.array([cal_partition_func(En, gn, T) for T in Ts], dtype=float)


def save_hitran_partition_func(hitran_df, Ntemp, Tmax):
    """
    Calculate and save partition functions for HITRAN database.

    Computes LTE partition functions at specified temperature intervals
    directly from HITRAN linelist data and saves results in .pf format file.

    Parameters
    ----------
    hitran_df : pd.DataFrame
        HITRAN DataFrame with linelist data
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
    Ts = np.array(range(Ntemp, Tmax+1, Ntemp)) 
    partition_func = calculate_linelist_partition_func(hitran_df, Ts)
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
