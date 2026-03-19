# Stick spectra and Cross sections

If you want to calculate stick spectra and cross sections simultaneously, you can set `StickSpectra` and `CrossSections` both to 1 and the program can calculate them together. This is more efficient than calculating them separately as it avoids redundant transition loading.

## Parameters

This calculation combines all parameters from both stick spectra and cross section calculations. For detailed descriptions of each parameter, please refer to [**Stick spectra**](`https://pyexocross.readthedocs.io/en/latest/stickspectra.html`) and [**Cross sections**](`https://pyexocross.readthedocs.io/en/latest/crosssections.html`).

1. **Temperatures and Wavenumber/Wavelength Ranges**: Shared parameters for both calculations.
2. **Pressures and Broadeners**: Applied to the cross section calculations.
3. **Line profiles**: Used for cross section profiles.
4. **Plotting**: Plots can be generated for both simultaneously using their respective configuration parameters (`PlotStickSpectra(Y/N)` and `PlotCrossSection(Y/N)`).

### Example Configuration

```bash
# Data source #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
SpeciesID                               666


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/MgH_ExoMol.log


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            1
CrossSections                           1


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                +/-  e/f   ElecState    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %12s         %3d   %3d      %5.1f    %5.1f


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            300,3000                  # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          +/-[]   e/f[]  ElecState[]  v[1,1;1,0;2,;,0]  


# Calculate cross sections #
Pressures                               1                         # Pressures in unit bar
Npoints/BinSize                         BinSize    0.1            # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt      
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule

# Plotting cross sections #
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30 cm2/molecule
```
