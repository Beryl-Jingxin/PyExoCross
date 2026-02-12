"""
Plotting functions for oscillator strengths.

This module provides functions for visualizing oscillator strength data.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ..base.utils import Timer, ensure_dir

# Plot Oscillator Strength
def plot_oscillator_strength(oscillator_strength_df):
    """
    Plot oscillator strength as a function of wavenumber or wavelength.

    Creates a plot of oscillator strength (gf or f) versus wavenumber
    or wavelength based on configuration settings.

    Parameters
    ----------
    oscillator_strength_df : pd.DataFrame
        DataFrame with columns 'v' (wavenumber) and 'os' (oscillator strength)
    """
    # Import legacy-style configuration variables from core (set via Config.to_globals()).
    from pyexocross.core import ( 
        save_path,
        data_info,
        database,
        PlotOscillatorStrengthMethod,
        PlotOscillatorStrengthWnWl,
        PlotOscillatorStrengthUnit,
        limitYaxisOS,
        gfORf,
    )

    print('Plotting oscillator strength ...')
    tp = Timer()
    tp.start()
    oscillator_strength_df = oscillator_strength_df[np.isfinite(oscillator_strength_df['os'])]
    oscillator_strength_df = oscillator_strength_df.reset_index(drop=True)
    v = oscillator_strength_df['v']
    osc_str = oscillator_strength_df['os']
    parameters = {'axes.labelsize': 18, 
                  'legend.fontsize': 18,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14}
    plt.rcParams.update(parameters)
    os_plot_folder = save_path + 'oscillator_strength/plots/'+data_info[0]+'/'+database+'/'
    ensure_dir(os_plot_folder)
    plt.figure(figsize=(12, 6))
    if PlotOscillatorStrengthMethod == 'LOG':
        plt.ylim([limitYaxisOS, 10*max(osc_str)])
        plt.semilogy()
    else:
        plt.ylim([limitYaxisOS, 1.05*max(osc_str)])
    #plt.title(database+' '+data_info[0]+' oscillator strengths') 
    if PlotOscillatorStrengthWnWl == 'WN':
        plt.plot(v, osc_str, label=data_info[0], linewidth=0.4)
        plt.xlabel('Wavenumber, cm⁻¹')
    elif PlotOscillatorStrengthWnWl == 'WL' and PlotOscillatorStrengthUnit == 'um':
        plt.plot(1e4/v, osc_str, label=data_info[0], linewidth=0.4)
        plt.xlabel('Wavelength, μm')
    elif PlotOscillatorStrengthWnWl == 'WL' and PlotOscillatorStrengthUnit == 'nm':
        plt.plot(1e7/v, osc_str, label=data_info[0], linewidth=0.4)
        plt.xlabel('Wavelength, nm') 
    else:
        raise ValueError('Please wirte the unit of wavelength in the input file: um or nm.')       
    plt.ylabel('Oscillator strength, ('+gfORf.lower()+')')
    plt.legend()
    leg = plt.legend()                  # Get the legend object.
    for line in leg.get_lines():
        line.set_linewidth(1.0)         # Change the line width for the legend.
    os_plot = (os_plot_folder+data_info[0]+'__'+data_info[1]+'__'+database+'__'+gfORf.lower()
               +'__'+PlotOscillatorStrengthWnWl.lower()+'__'+PlotOscillatorStrengthUnit.lower()+'__os.png')
    plt.savefig(os_plot, dpi=500)
    plt.show()
    tp.end()
    print('Oscillator strengths plot has been saved:', os_plot)  
    