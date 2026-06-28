"""
Plotting functions for PyExoCross.

This module provides functions for visualizing:
- Stick spectra
- Cross sections
- Oscillator strengths
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from ..base.utils import Timer, ensure_dir
from ..process.stick_xsec_filepath import temperature_string_base
from .mpl_safety import configure_mpl_large_path_rendering


def _axis_values_from_wavenumber(v, values, axis_kind, axis_unit):
    """Return paired x/value arrays in the requested axis unit, sorted by x."""
    v = np.asarray(v)
    values = np.asarray(values)
    if axis_kind == 'WN':
        x_value = v
        value_out = values
        unit_fn = 'cm-1__'
        axis_label = 'Wavenumber, cm⁻¹'
    elif axis_kind == 'WL' and axis_unit == 'um':
        if np.any(v <= 0):
            print(
                'Warning: Wavenumber values at or below 0 cm⁻¹ cannot be '
                'converted to a finite wavelength and were omitted. Set '
                'min_range above 0 when plotting in wavelength.'
            )
        mask = v > 0
        x_value = 1e4 / v[mask]
        value_out = values[mask]
        unit_fn = 'um__'
        axis_label = 'Wavelength, μm'
    elif axis_kind == 'WL' and axis_unit == 'nm':
        if np.any(v <= 0):
            print(
                'Warning: Wavenumber values at or below 0 cm⁻¹ cannot be '
                'converted to a finite wavelength and were omitted. Set '
                'min_range above 0 when plotting in wavelength.'
            )
        mask = v > 0
        x_value = 1e7 / v[mask]
        value_out = values[mask]
        unit_fn = 'nm__'
        axis_label = 'Wavelength, nm'
    else:
        raise ValueError('Please wirte the unit of wavelength in the input file: um or nm.')

    valid = np.isfinite(x_value) & np.isfinite(value_out)
    x_value = x_value[valid]
    value_out = value_out[valid]
    order = np.argsort(x_value)
    return x_value[order], value_out[order], unit_fn, axis_label


# Plot stick spectra
def plot_stick_spectra(stick_spectra_df, T=None, Tvib=None, Trot=None):
    """
    Plot stick spectra (line intensities) as vertical lines.

    Creates a plot of stick spectra with vertical lines representing
    transition intensities at their wavenumber positions.

    Parameters
    ----------
    stick_spectra_df : pd.DataFrame
        DataFrame with columns 'v' (wavenumber) and 'S' (intensity)
    T : float, optional
        Temperature in Kelvin for LTE or population method
    Tvib : float, optional
        Vibrational temperature in Kelvin for non-LTE two-temperature method
    Trot : float, optional
        Rotational temperature in Kelvin for non-LTE methods
    """  
    from pyexocross.core import ( 
        save_path,
        data_info,
        database,
        PlotStickSpectraWnWl,
        PlotStickSpectraUnit,
        PlotStickSpectraMethod,
        limitYaxisStickSpectra,
        wn_wl,
        wn_wl_unit,
        abs_emi,
        LTE_NLTE,
        photo,
        NLTEMethod,
        UncFilter,
        threshold,
        T_list,
        Tvib_list,
        Trot_list,
    )
    print('\nPlotting stick spectra ...')
    tp = Timer()
    tp.start()
    stick_spectra_df_plot = stick_spectra_df[np.isfinite(stick_spectra_df['S'])]
    stick_spectra_df_plot = stick_spectra_df_plot.reset_index(drop=True)
    configure_mpl_large_path_rendering()
    
    if len(stick_spectra_df_plot) == 0:
        print('Warning: No valid stick spectra data to plot. Skipping plot.')
        return
    
    v = stick_spectra_df_plot['v'].values
    S = stick_spectra_df_plot['S'].values
    del stick_spectra_df_plot
    min_v = min(v)
    max_v = max(v)
    
    parameters = {'axes.labelsize': 18, 
                  'legend.fontsize': 18,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14}
    plt.rcParams.update(parameters)
    ss_plot_folder = save_path + 'stick_spectra/plots/'+data_info[0]+'/'+database+'/'
    ensure_dir(ss_plot_folder)
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Store original S for y-axis limits (before any filtering)
    S_original = S.copy()
    v_value, S, unit_fn, axis_label = _axis_values_from_wavenumber(
        v, S, PlotStickSpectraWnWl, PlotStickSpectraUnit
    )
    ax.set_xlabel(axis_label)
    if len(v_value) > 0:
        ax.set_xlim([np.min(v_value), np.max(v_value)])
    else:
        print(f'Warning: No valid data after filtering. Skipping plot.')
        plt.close()
        return
    # Check if we have data after filtering
    if len(v_value) == 0 or len(S) == 0:
        print(f'Warning: No valid data to plot after filtering (v length: {len(v_value)} {wn_wl_unit}, S length: {len(S)}) . Skipping plot.')
        plt.close()
        return
    
    if abs_emi == 'Ab':
        stick_y_unit = 'cm/molecule'
    else:
        stick_y_unit = 'erg (s molecule sr)⁻¹'

    # Debug: Check data ranges
    if len(v_value) > 0 and len(S) > 0:
        plot_unit = PlotStickSpectraUnit.replace('um', 'μm').replace('cm-1', 'cm⁻¹')
        print(f'Info: Plotting stick spectra {len(v_value)} points, x range: [{np.min(v_value):.2f}, {np.max(v_value):.2f}] {plot_unit}, S range: [{np.min(S):.2e}, {np.max(S):.2e}] {stick_y_unit}')
    
    # Recalculate y-axis limits based on filtered S if needed
    if PlotStickSpectraMethod == 'LOG':
        max_S_val = max(S) if len(S) > 0 else max(S_original)
        ax.set_ylim([limitYaxisStickSpectra, 10*max_S_val])
        ax.semilogy()
    else:
        max_S_val = max(S) if len(S) > 0 else max(S_original)
        ax.set_ylim([limitYaxisStickSpectra, 1.05*max_S_val])
    
    ax.set_ylabel(f'Intensity, {stick_y_unit}')
    # T, Tvib, Trot should be passed as parameters
    T_val = None
    Tvib_val = None
    Trot_val = None
    if NLTEMethod == 'L' or NLTEMethod == 'P':
        T_val = T if T is not None else (T_list[0] if len(T_list) > 0 else T_list)
        ax.vlines(v_value, 0, S, label=f'T = {T_val} K', linewidth=1.5, alpha=1)
    elif NLTEMethod == 'T':
        T_val = T if T is not None else (T_list[0] if len(T_list) > 0 else T_list)
        Tvib_val = Tvib if Tvib is not None else (Tvib_list[0] if len(Tvib_list) > 0 else Tvib_list)
        Trot_val = Trot if Trot is not None else (Trot_list[0] if len(Trot_list) > 0 else Trot_list)
        ax.vlines(v_value, 0, S, label=f'T$_{{vib}}$ = {Tvib_val} K, T$_{{rot}}$ = {Trot_val} K', linewidth=1.5, alpha=1)
    elif NLTEMethod == 'D':
        Trot_val = Trot if Trot is not None else (Trot_list[0] if len(Trot_list) > 0 else Trot_list)
        ax.vlines(v_value, 0, S, label=f'T$_{{rot}}$ = {Trot_val} K', linewidth=1.5, alpha=1)
    else:
        raise ValueError("Please choose 'LTE' or 'Non-LTE'; if choose 'Non-LTE', please choose one non-LTE method from: 'T', 'D' or 'P'.")
    #plt.title(database+' '+data_info[0]+' intensity') 
    leg = ax.legend()                  # Get the legend object.
    # for line in leg.legend_handles:
    #     line.set_height(1.5)           # Change the line width for the legend.
    plot_min_v = np.min(v_value)
    plot_max_v = np.max(v_value)
    str_min_wnl = str(int(np.floor(plot_min_v)))
    str_max_wnl = str(int(np.ceil(plot_max_v)))
    T_str = temperature_string_base(T_val, Tvib_val, Trot_val, NLTEMethod)
    ss_plot = (ss_plot_folder + '__'.join(data_info) + '__' + T_str + '__' + PlotStickSpectraWnWl.lower()
               + str_min_wnl + '-' + str_max_wnl + unit_fn
               + 'unc' + str(UncFilter) + '__thres' + str(threshold)
               + '__' + database + '__' + abs_emi + photo + LTE_NLTE + '.png')
    plt.savefig(ss_plot, dpi=500)
    plt.close()  # Close figure to free memory and ensure it's saved
    tp.end()
    print('Stick spectra plot has been saved:', ss_plot, '\n')
