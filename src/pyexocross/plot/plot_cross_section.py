"""
Plotting functions for cross sections.

This module provides functions for plotting and saving cross-section data.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ..base.utils import Timer, ensure_dir
from ..process.stick_xsec_filepath import (
    cross_section_filepath,
    crosssectiondetails,
    temperature_pressure_string,
)
from .mpl_safety import configure_mpl_large_path_rendering


def _axis_values_from_wavenumber(wn, values, axis_kind, axis_unit):
    """Return paired x/value arrays in the requested axis unit, sorted by x."""
    wn = np.asarray(wn)
    values = np.asarray(values)
    if axis_kind == 'WN':
        x_value = wn
        value_out = values
        unit_fn = 'cm-1__'
        axis_label = 'Wavenumber, cm⁻¹'
    elif axis_kind == 'WL' and axis_unit == 'um':
        if np.any(wn <= 0):
            print(
                'Warning: Wavenumber values at or below 0 cm⁻¹ cannot be '
                'converted to a finite wavelength and were omitted. Set '
                'min_range above 0 when plotting in wavelength.'
            )
        mask = wn > 0
        x_value = 1e4 / wn[mask]
        value_out = values[mask]
        unit_fn = 'um__'
        axis_label = 'Wavelength, μm'
    elif axis_kind == 'WL' and axis_unit == 'nm':
        if np.any(wn <= 0):
            print(
                'Warning: Wavenumber values at or below 0 cm⁻¹ cannot be '
                'converted to a finite wavelength and were omitted. Set '
                'min_range above 0 when plotting in wavelength.'
            )
        mask = wn > 0
        x_value = 1e7 / wn[mask]
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


# Plot and Save Results
def save_xsec_file_plot(wn, xsec, database, profile_label, T=None, P=None, temp_idx=None,
                        Tvib_list_param=None, Trot_list_param=None):
    """
    Save cross-section data to file and optionally create plots.

    Saves cross-section data in .xsec format and optionally generates
    plots in wavenumber or wavelength units based on configuration.

    Parameters
    ----------
    wn : np.ndarray
        Wavenumber array, shape (n_points,)
    xsec : np.ndarray
        Cross-section array, shape (n_points,)
    database : str
        Database name (e.g., 'ExoMol', 'HITRAN')
    profile_label : str
        Line profile label (e.g., 'Voigt', 'Doppler')
    T : float, optional
        Temperature in Kelvin (for single temperature)
    P : float or list, optional
        Pressure in bar (for single or range)
    temp_idx : int, optional
        Temperature index for non-LTE two-temperature method
    """             
    from pyexocross.core import (  # type: ignore
        save_path,
        data_info,
        NLTEMethod,
        wn_wl,
        wn_wl_unit,
        UncFilter,
        threshold,
        abs_emi,
        PlotCrossSectionYN,
        PlotCrossSectionWnWl,
        PlotCrossSectionUnit,
        PlotCrossSectionMethod,
        limitYaxisXsec,
        LTE_NLTE,
        photo,
        bin_size,
        cutoff,
        T_list,
    )
    xsecs_foldername = save_path+'xsecs/files/'+data_info[0]+'/'+database+'/'
    ensure_dir(xsecs_foldername)
    configure_mpl_large_path_rendering()
    wn = np.asarray(wn)
    xsec = np.asarray(xsec)
    min_v = min(wn)
    max_v = max(wn)
    # Check if profile is pressure-dependent (Lorentzian/Voigt need pressure, but Doppler/Gaussian don't)
    pressure_dependent = profile_label not in ['Gaussian', 'Doppler', 'Binned Doppler', 'Binned Gaussion']
    if wn_wl == 'WN':
        bin_unit_fn = 'cm-1__'
    elif wn_wl == 'WL' and wn_wl_unit == 'um':
        bin_unit_fn = 'um__'
    elif wn_wl == 'WL' and wn_wl_unit == 'nm':
        bin_unit_fn = 'nm__'
    else:
        raise ValueError('Please choose wavenumber or wavelength and type in correct format: wn or wl.')

    # Resolve local Tvib/Trot lists (allow explicit args or fall back to globals)
    if Tvib_list_param is not None:
        Tvib_list_use = Tvib_list_param
    else:
        Tvib_list_use = globals().get('Tvib_list', [])
    if Trot_list_param is not None:
        Trot_list_use = Trot_list_param
    else:
        Trot_list_use = globals().get('Trot_list', [])

    # Build temperature (and optional pressure) string using shared helper
    T_str = temperature_pressure_string(
        T=T,
        P=P,
        temp_idx=temp_idx,
        NLTEMethod=NLTEMethod,
        Tvib_list=Tvib_list_use,
        Trot_list=Trot_list_use,
        pressure_dependent=pressure_dependent,
    )
    if 'L' not in wn_wl:
        # Save cross sections into .xsec file.
        print('\nSaving cross sections into file ...')   
        ts = Timer()    
        ts.start()
        xsec_df = pd.DataFrame()
        xsec_df['wavenumber'] = wn
        xsec_df['cross-section'] = xsec
        str_min_wn = str(int(np.floor(min_v)))
        str_max_wn = str(int(np.ceil(max_v)))
        xsec_filepath = cross_section_filepath(
            xsecs_foldername,
            data_info,
            T,
            P,
            temp_idx,
            Tvib_list_use,
            Trot_list_use,
            str_min_wn,
            str_max_wn,
            'cm-1__',
            wn_wl,
            UncFilter,
            threshold,
            database,
            abs_emi,
            bin_size,
            cutoff,
            profile_label,
            LTE_NLTE,
            photo,
            NLTEMethod,
            pressure_dependent,
        )
        np.savetxt(xsec_filepath, xsec_df, fmt="%15.6f %15.8E")
        ts.end()
        # File info printed once by caller
        print('Cross sections file has been saved:', xsec_filepath)
        
        if abs_emi == 'Ab':
            xsec_y_unit = 'cm²/molecule'
        else:
            xsec_y_unit = 'erg cm (s molecule sr)⁻¹'

        if PlotCrossSectionYN == 'Y':
            print('\nPlotting cross sections ...')
            tp = Timer()
            tp.start()
            plots_foldername = save_path+'xsecs/plots/'+data_info[0]+'/'+database+'/'
            ensure_dir(plots_foldername)
            #plt.legend(fancybox=True, framealpha=0.0)
            parameters = {'axes.labelsize': 18,
                          'legend.fontsize': 18,
                          'xtick.labelsize': 14,
                          'ytick.labelsize': 14}
            plt.rcParams.update(parameters)
            # Plot cross sections and save it as .png.
            plt.figure(figsize=(12, 6))
            v_value, xsec_plot, unit_fn, axis_label = _axis_values_from_wavenumber(
                wn, xsec, PlotCrossSectionWnWl, PlotCrossSectionUnit
            )
            plt.xlabel(axis_label)
            if len(v_value) > 0:
                plot_min_v = np.min(v_value)
                plot_max_v = np.max(v_value)
                plt.xlim([plot_min_v, plot_max_v])
            else:
                plot_min_v = min_v
                plot_max_v = max_v
            
            # Check if we have data to plot
            if len(xsec_plot) == 0 or len(v_value) == 0:
                print(f'Warning: No valid data to plot for cross sections (xsec_plot length: {len(xsec_plot)}, v_value length: {len(v_value)}). Original xsec length: {len(xsec)}, non-zero xsec count: {np.sum(xsec > 0) if len(xsec) > 0 else 0}. Skipping plot.')
                plt.close()
                return
            
            # Check if xsec_plot has any non-zero values
            if len(xsec_plot) > 0:
                if np.all(~np.isfinite(xsec_plot)):
                    print(f'Warning: All cross section values are invalid (NaN/Inf). Skipping plot.')
                    plt.close()
                    return
                elif np.all(xsec_plot == 0):
                    print(f'Warning: All cross section values are zero. This might indicate no transitions in the range. Skipping plot.')
                    plt.close()
                    return
                else:
                    # Debug: Check data ranges
                    plot_unit = PlotCrossSectionUnit.replace('um', 'μm').replace('cm-1', 'cm⁻¹')
                    print(f'Info: Plotting cross sections: {len(v_value)} points, x range: [{np.min(v_value):.2f}, {np.max(v_value):.2f}] {plot_unit}, xsec range: [{np.min(xsec_plot):.2e}, {np.max(xsec_plot):.2e}] {xsec_y_unit}, non-zero xsec count: {np.sum(xsec_plot > 0)}')
            else:
                print(f'Warning: xsec_plot is empty after filtering. Original xsec length: {len(xsec)}, non-zero xsec count: {np.sum(xsec > 0) if len(xsec) > 0 else 0}. Skipping plot.')
                plt.close()
                return
            
            # Ensure we have valid data before plotting
            if len(v_value) != len(xsec_plot):
                print(f'Error: Mismatch in array lengths: v length={len(v_value)}, xsec length={len(xsec_plot)}. Skipping plot.')
                plt.close()
                return
            
            # Add pressure to label if profile is pressure-dependent
            P_label = ''
            if pressure_dependent and P is not None:
                P_label = f', P = {P} bar'
            
            if NLTEMethod == 'L' or NLTEMethod == 'P':
                label_str = f'T = {T} K{P_label}, {profile_label}' if T is not None else f'T = {min(T_list)}-{max(T_list)} K{P_label}, {profile_label}'
                plt.plot(v_value, xsec_plot, label=label_str, linewidth=0.4)  
            elif NLTEMethod == 'T':
                if T is not None and temp_idx is not None:
                    Tvib_val = Tvib_list_use[temp_idx]
                    Trot_val = Trot_list_use[temp_idx]
                    label_str = f'T$_{{vib}}$ = {Tvib_val} K, T$_{{rot}}$ = {Trot_val} K{P_label}, {profile_label}'
                else:
                    label_str = (
                        f'T$_{{vib}}$ = {min(Tvib_list_use)}-{max(Tvib_list_use)} K, '
                        f'T$_{{rot}}$ = {min(Trot_list_use)}-{max(Trot_list_use)} K{P_label}, {profile_label}'
                    )
                plt.plot(v_value, xsec_plot, label=label_str, linewidth=0.4)
            elif NLTEMethod == 'D':
                if T is not None and temp_idx is not None:
                    Trot_val = Trot_list_use[temp_idx]
                    label_str = f'T$_{{rot}}$ = {Trot_val} K{P_label}, {profile_label}'
                else:
                    label_str = (
                        f'T$_{{rot}}$ = {min(Trot_list_use)}-{max(Trot_list_use)} K{P_label}, {profile_label}'
                    )
                plt.plot(v_value, xsec_plot, label=label_str, linewidth=0.4)
            else:
                raise ValueError("Please choose 'LTE' or 'Non-LTE'; if choose 'Non-LTE', please choose one non-LTE method from: 'T', 'D' or 'P'.")
            if PlotCrossSectionMethod == 'LOG':
                # Use nanmax to ignore NaN values
                max_xsec_val = np.nanmax(xsec_plot) if len(xsec_plot) > 0 else np.nanmax(xsec)
                # Check if max_xsec_val is valid (not NaN or Inf)
                if np.isfinite(max_xsec_val) and max_xsec_val > 0:
                    plt.ylim([limitYaxisXsec, 10*max_xsec_val])
                else:
                    print(f"Warning: Invalid max_xsec_val ({max_xsec_val}), skipping ylim setting")
                plt.semilogy()
            else:
                # Use nanmax to ignore NaN values
                max_xsec_val = np.nanmax(xsec_plot) if len(xsec_plot) > 0 else np.nanmax(xsec)
                # Check if max_xsec_val is valid (not NaN or Inf)
                if np.isfinite(max_xsec_val) and max_xsec_val > 0:
                    plt.ylim([limitYaxisXsec, 1.05*max_xsec_val])
                else:
                    print(f"Warning: Invalid max_xsec_val ({max_xsec_val}), skipping ylim setting")
            # plt.title(database+' '+data_info[0]+' '+abs_emi+' Cross-Section with '+ profile_label)
            plt.ylabel(f'Cross-section, {xsec_y_unit}')
            plt.legend()
            leg = plt.legend()                  # Get the legend object.
            for line in leg.get_lines():
                line.set_linewidth(1.0)         # Change the line width for the legend.
            # Use appropriate values for filename based on plot type
            str_min_wnl = str(int(np.floor(plot_min_v)))
            str_max_wnl = str(int(np.ceil(plot_max_v)))
            xsec_plotpath = (plots_foldername + '__'.join(data_info) + '__' + T_str + '__' + PlotCrossSectionWnWl.lower()
                             + str_min_wnl + '-' + str_max_wnl + unit_fn + 'unc' + str(UncFilter)
                             + '__thres' + str(threshold) + '__' + database + '__' + abs_emi
                             + crosssectiondetails(bin_size, bin_unit_fn, cutoff, profile_label)
                             + photo + LTE_NLTE + '.png')
            plt.savefig(xsec_plotpath, dpi=500)
            plt.close()  # Close figure to free memory and ensure it's saved
            tp.end()
            print('Cross sections plot has been saved:', xsec_plotpath)
        else:
            pass
    elif 'L' in wn_wl:
        if wn_wl_unit == 'um':
            unit_change = 1e4
            unit_str = 'μm'
            unit_ffn = 'um__'
        elif wn_wl_unit == 'nm':
            unit_change = 1e7
            unit_str = 'nm'
            unit_ffn = 'nm__'
        else:
            raise ValueError('Please wirte the unit of wavelength in the input file: um or nm.')   
        wl, xsec_save, _, _ = _axis_values_from_wavenumber(wn, xsec, 'WL', wn_wl_unit)
        min_wl = min(wl)
        max_wl = max(wl)
        # min_wl = '%.02f' % min(wl)
        # max_wl = '%.02f' % max(wl)
        # Save cross sections into .xsec file.
        print('\nSaving cross sections into file ...')
        ts = Timer()    
        ts.start()
        xsec_df = pd.DataFrame()
        xsec_df['wavelength'] = wl
        xsec_df['cross-section'] = xsec_save
        str_min_wl = str(int(np.floor(min_wl)))
        str_max_wl = str(int(np.ceil(max_wl)))
        xsec_filepath = cross_section_filepath(
            xsecs_foldername,
            data_info,
            T,
            P,
            temp_idx,
            Tvib_list_use,
            Trot_list_use,
            str_min_wl,
            str_max_wl,
            unit_ffn,
            wn_wl,
            UncFilter,
            threshold,
            database,
            abs_emi,
            bin_size,
            cutoff,
            profile_label,
            LTE_NLTE,
            photo,
            NLTEMethod,
            pressure_dependent,
        )
        np.savetxt(xsec_filepath, xsec_df, fmt="%15.8E %15.8E")
        ts.end()
        # File info printed once by caller
        print('Cross sections file has been saved:', xsec_filepath)       
        
        if abs_emi == 'Ab':
            xsec_y_unit = 'cm²/molecule'
        else:
            xsec_y_unit = 'erg cm (s molecule sr)⁻¹'

        if PlotCrossSectionYN == 'Y':
            print('\nPlotting cross sections ...')
            tp = Timer()
            tp.start()
            plots_foldername = save_path+'xsecs/plots/'+data_info[0]+'/'+database+'/'
            ensure_dir(plots_foldername)
            # plt.legend(fancybox=True, framealpha=0.0)
            parameters = {'axes.labelsize': 18, 
                          'legend.fontsize': 18,
                          'xtick.labelsize': 14,
                          'ytick.labelsize': 14}
            plt.rcParams.update(parameters)
            # Plot cross sections and save it as .png.
            plt.figure(figsize=(12, 6))
            v_value, xsec_plot, unit_pfn, plot_unit_str = _axis_values_from_wavenumber(
                wn, xsec, PlotCrossSectionWnWl, PlotCrossSectionUnit
            )
            if len(v_value) == 0:
                print('Warning: No valid data to plot for cross sections. Skipping plot.')
                plt.close()
                return
            min_v_value = np.min(v_value)
            max_v_value = np.max(v_value)
            plt.xlim([min_v_value, max_v_value])

            P_label = ''
            if pressure_dependent and P is not None:
                P_label = f', P = {P} bar'

            if NLTEMethod == 'L' or NLTEMethod == 'P':
                label_str = f'T = {T} K{P_label}, {profile_label}' if T is not None else f'T = {min(T_list)}-{max(T_list)} K{P_label}, {profile_label}'
                plt.plot(v_value, xsec_plot, label=label_str, linewidth=0.4)
            elif NLTEMethod == 'T':
                if T is not None and temp_idx is not None:
                    Tvib_val = Tvib_list_use[temp_idx]
                    Trot_val = Trot_list_use[temp_idx]
                    label_str = f'T$_{{vib}}$ = {Tvib_val} K, T$_{{rot}}$ = {Trot_val} K{P_label}, {profile_label}'
                else:
                    label_str = (
                        f'T$_{{vib}}$ = {min(Tvib_list_use)}-{max(Tvib_list_use)} K, '
                        f'T$_{{rot}}$ = {min(Trot_list_use)}-{max(Trot_list_use)} K{P_label}, {profile_label}'
                    )
                plt.plot(v_value, xsec_plot, label=label_str, linewidth=0.4)
            elif NLTEMethod == 'D':
                if T is not None and temp_idx is not None:
                    Trot_val = Trot_list_use[temp_idx]
                    label_str = f'T$_{{rot}}$ = {Trot_val} K{P_label}, {profile_label}'
                else:
                    label_str = (
                        f'T$_{{rot}}$ = {min(Trot_list_use)}-{max(Trot_list_use)} K{P_label}, {profile_label}'
                    )
                plt.plot(v_value, xsec_plot, label=label_str, linewidth=0.4)
            else:
                raise ValueError("Please choose 'LTE' or 'Non-LTE'; if choose 'Non-LTE', please choose one non-LTE method from: 'T', 'D' or 'P'.")          
            if PlotCrossSectionMethod == 'LOG':
                plt.ylim([limitYaxisXsec, 10*max(xsec_plot)])
                plt.semilogy()
            else:
                plt.ylim([limitYaxisXsec, 1.05*max(xsec_plot)])
            #plt.title(database+' '+data_info[0]+' '+abs_emi+' Cross-Section with '+ profile_label) 
            plt.xlabel(plot_unit_str)
            plt.ylabel(f'Cross-section, {xsec_y_unit}')
            plt.legend()
            leg = plt.legend()                  # Get the legend object.
            for line in leg.get_lines():
                line.set_linewidth(1.0)         # Change the line width for the legend.
            str_min_wnl = str(int(np.floor(min_v_value)))
            str_max_wnl = str(int(np.ceil(max_v_value)))
            xsec_plotpath = (plots_foldername+'__'.join(data_info)+'__'+T_str+'__'+PlotCrossSectionWnWl.lower()
                             +str_min_wnl+'-'+str_max_wnl+unit_pfn+'unc'+str(UncFilter)
                             +'__thres'+str(threshold)+'__'+database+'__'+abs_emi 
                             + crosssectiondetails(bin_size, unit_ffn, cutoff, profile_label)
                             + photo + LTE_NLTE + '.png')
            plt.savefig(xsec_plotpath, dpi=500)
            plt.close()  # Close figure to free memory and ensure it's saved
            tp.end()
            print('Cross sections plot has been saved:', xsec_plotpath)
        else:
            pass
    else:
        raise ValueError('Please choose wavenumber or wavelength and type in correct format: wn or wl.')
