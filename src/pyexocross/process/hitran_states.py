"""
Build approximate state lists from HITRAN transition records.

HITRAN ``.par`` files are transition lists rather than complete state files.
The helpers here collapse repeated transition records into the best available
state list for calculations that only need state energies and degeneracies.
"""
import pandas as pd


def normalize_state_symmetry(value):
    """
    Normalize HITRAN lower-state symmetry labels for state grouping.

    Some HITRAN class-1 local QN fields store two parity/symmetry characters
    such as "ef" or "fe". For a lower-state key, the second character is the
    lower-state symmetry; single-character labels are already state labels.
    """
    text = str(value).strip()
    if len(text) == 2 and all(char in 'ef+-' for char in text):
        return text[1]
    return text


def build_hitran_partition_states(hitran_linelist_df):
    """
    Build approximate HITRAN states for state-sum calculations.

    Hyperfine components and multiple branches can make one physical lower
    state appear many times. Collapse records by lower-state quantum numbers
    excluding F, use the lowest listed component energy, and sum the unique F
    statistical weights for that physical state.
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


def build_hitran_states(hitran_df):
    from pyexocross.process.hitran_qn import hitran_linelist_QN

    hitran_linelist_df, _ = hitran_linelist_QN(hitran_df)
    return build_hitran_partition_states(hitran_linelist_df)


def hitran_state_arrays(hitran_df):
    states_df = build_hitran_states(hitran_df)
    En = states_df['E'].astype('float').values
    gn = states_df['g'].astype('int').values
    return En, gn
