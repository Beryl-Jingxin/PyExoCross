"""
Filter linelist by quantum numbers.

This module provides functions for filtering molecular line lists based on
quantum number values.
"""
import pandas as pd
from itertools import chain

def QNfilter_linelist(linelist_df, QNs_value, QNs_label):
    """
    Filter linelist DataFrame based on quantum number values.

    Filters transitions by matching quantum number values for upper and lower
    states. Handles wildcard matching where empty strings in QNs_value match
    all values for that quantum number.

    Parameters
    ----------
    linelist_df : pd.DataFrame
        Linelist DataFrame with quantum number columns (with ' and " suffixes)
    QNs_value : list of list of str
        List of quantum number value lists, each containing values to match
    QNs_label : list of str
        List of quantum number labels to filter on

    Returns
    -------
    pd.DataFrame
        Filtered linelist DataFrame
    """
    for i in range(len(QNs_label)):
        if QNs_value[i] != ['']:
            linelist_df[i] = linelist_df[[QNs_label[i]+"'",QNs_label[i]+'"']].fillna('').agg(','.join, axis=1).astype(str).str.replace(' ', '')
            uval = linelist_df[QNs_label[i]+"'"].astype(str).str.replace(' ', '').drop_duplicates().values
            lval = linelist_df[QNs_label[i]+'"'].astype(str).str.replace(' ', '').drop_duplicates().values
            vallist = []
            ulist = []
            llist = []
            for qnval in QNs_value[i]:
                if '' not in qnval.split(','):
                    vallist.append(qnval)
                elif qnval.split(',')[0] == '':
                    ulist.append([qnval.replace(",",val+",") for val in uval])
                elif qnval.split(',')[1] == '':
                    llist.append([qnval.replace(",",","+val) for val in lval])
            QNs_value[i] = vallist+list(chain(*ulist))+list(chain(*llist))
            linelist_df = linelist_df[linelist_df[i].isin(QNs_value[i])].drop(columns=[i])   
    # linelist_df = linelist_df.sort_values('v')  
    return(linelist_df)   
