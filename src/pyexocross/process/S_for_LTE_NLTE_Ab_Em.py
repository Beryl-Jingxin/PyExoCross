from pyexocross.calculation.calculate_intensity import (
    cal_abscoefs,
    cal_abscoefs_nlte_2T,
    cal_abscoefs_nlte_nvib,
    cal_abscoefs_nlte_pop,
)
from pyexocross.calculation.calculate_emissivity import (
    cal_emicoefs,
    cal_emicoefs_nlte_2T,
    cal_emicoefs_nlte_nvib,
    cal_emicoefs_nlte_pop,
)
from pyexocross.process.T_n_val import get_temp_vals
def S_for_LTE_NLTE_Ab_Em(ll_df,T_list,Tvib_list,Trot_list,Q_arr,abs_emi,NLTEMethod,temp_idx=None):
    """
    Calculate the stick spectra for LTE and NLTE.
    """
    from pyexocross.core import abundance
    A = ll_df['A'].values
    v = ll_df['v'].values
    T, Tvib, Trot = get_temp_vals(temp_idx, NLTEMethod, T_list, Tvib_list, Trot_list)
    if abs_emi == 'Ab':
        if NLTEMethod == 'L':
            # Process single temperature if temp_idx is provided
            if temp_idx is not None:
                S_arr = cal_abscoefs([T],[Q_arr[temp_idx]],
                                     ll_df['E"'].values,
                                     ll_df["g'"].values,
                                     A,v,abundance)
            else:
                S_arr = cal_abscoefs(T_list,Q_arr,
                                     ll_df['E"'].values,
                                     ll_df["g'"].values,
                                     A,v,abundance)
        elif NLTEMethod == 'T':
            if temp_idx is not None:
                S_arr = cal_abscoefs_nlte_2T([Tvib],[Trot],[Q_arr[temp_idx]],
                                             ll_df['Evib"'],
                                             ll_df['Erot"'],
                                             ll_df["g'"],
                                             A,v,abundance)
            else:
                S_arr = cal_abscoefs_nlte_2T(Tvib_list,Trot_list,Q_arr,       
                                             ll_df['Evib"'],
                                             ll_df['Erot"'],
                                             ll_df["g'"],
                                             A,v,abundance)    
        elif NLTEMethod == 'D':
            if temp_idx is not None:
                S_arr = cal_abscoefs_nlte_nvib([T],[Trot],[Q_arr[temp_idx]],
                                                ll_df['nvib"'],
                                                ll_df['E"'],
                                                ll_df["g'"],
                                                A,v,abundance)
            else:
                S_arr = cal_abscoefs_nlte_nvib(T_list,Trot_list,Q_arr,
                                               ll_df['nvib"'],
                                               ll_df['E"'],
                                               ll_df["g'"],
                                               A,v,abundance)
        elif NLTEMethod == 'P':
            if temp_idx is not None:
                S_arr = cal_abscoefs_nlte_pop([T],
                                               ll_df['pop"'],
                                               ll_df["g'"],
                                               ll_df['g"'],
                                               A,v,abundance)
            else:
                S_arr = cal_abscoefs_nlte_pop(T_list,
                                              ll_df['pop"'],
                                              ll_df["g'"],
                                              ll_df['g"'],
                                              A,v,abundance)
        else:
            raise ValueError("Please choose 'LTE' or 'Non-LTE'; if choose 'Non-LTE', please choose one non-LTE method from: 'T', 'D' or 'P'.")
    elif abs_emi == 'Em':   
        if NLTEMethod == 'L':
            if temp_idx is not None:
                S_arr = cal_emicoefs([T],[Q_arr[temp_idx]],
                                     ll_df["E'"],
                                     ll_df["g'"],
                                     A,v,abundance)
            else:
                S_arr = cal_emicoefs(T_list,Q_arr,
                                     ll_df["E'"],
                                     ll_df["g'"],
                                     A,v,abundance)
        elif NLTEMethod == 'T':
            if temp_idx is not None:
                S_arr = cal_emicoefs_nlte_2T([Tvib],[Trot],[Q_arr[temp_idx]],
                                             ll_df["Evib'"],
                                             ll_df["Erot'"],
                                             ll_df["g'"],
                                             A,v,abundance)
            else:
                S_arr = cal_emicoefs_nlte_2T(Tvib_list,Trot_list,Q_arr,
                                             ll_df["Evib'"],
                                             ll_df["Erot'"],
                                             ll_df["g'"],
                                             A,v,abundance)
        elif NLTEMethod == 'D':   
            if temp_idx is not None:
                S_arr = cal_emicoefs_nlte_nvib([Trot],[Q_arr[temp_idx]],
                                               ll_df["nvib'"],
                                               ll_df["E'"],
                                               ll_df["g'"],
                                               A,v,abundance)
            else:
                S_arr = cal_emicoefs_nlte_nvib(Trot_list,Q_arr,
                                               ll_df["nvib'"],
                                               ll_df["E'"],
                                               ll_df["g'"],
                                               A,v,abundance)
        elif NLTEMethod == 'P':
            # Population method doesn't depend on temperature, so no temp_idx needed
            S_arr = cal_emicoefs_nlte_pop(ll_df["pop'"],
                                          A,v,abundance)
        else:
            raise ValueError("Please choose 'LTE' or 'Non-LTE'; if choose 'Non-LTE', please choose one non-LTE method from: 'T', 'D' or 'P'.")
    else:
        raise ValueError("Please choose one from: 'Absorption' or 'Emission'.") 
    return S_arr[0, :]
