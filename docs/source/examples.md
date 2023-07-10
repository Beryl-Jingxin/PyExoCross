# Examples

The examples of the whole input files for the ExoMol and HITRAN databases.

## Example for the ExoMol database

```bash
# Data source #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
MolIsoID                                501


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
StickSpectra                            0
CrossSections                           1


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f


# Conversion #
ConversionFormat                        1  
ConversionFrequncyRange                 0                 30000      
GlobalQNLabel                           eS       v        Omega
GlobalQNFormat                          %9s     %2d      %4s
LocalQNLabel                            J        e/f
LocalQNFormat                           %5.1f    %2s
ConvUncFilter(Y/N)                      N          0.01           # If Y, default value 0.001
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30
                           

# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
None


# Calculate stick spectra or cross sections #
Temperature                             300
Range                                   0          30000
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.001          # If Y, default value 0.001
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   N


# Calculate cross sections #
Pressure                                1
Npoints/BinSize                         Npoints    10001
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt        
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
Cutoff(Y/N)                             Y          100            # If Y, default value 25
QNsFilter(Y/N)                          N          par[+]   e/f[e]   v[0,1,2,3]  
DopplerHWHM(Y/N)                        Y          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   N
```

## Example for the HITRAN database

```bash
# Basic information #
Database                                HITRAN
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
MolIsoID                                501


# File path #
ReadPath                                /home/jingxin/data/pyexocross/conversion/24Mg-1H__XAB.par
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              1
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
StickSpectra                            0
CrossSections                           0


# Quantum numbers #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %1d      %7s      %7s


# Conversion #
ConversionFormat                        2  
ConversionUncertainty                   0.005
ConversionFrequncyRange                 0                 30000      
GlobalQNLabel                           eS       v        Omega
GlobalQNFormat                          %10s     %1d      %4s
LocalQNLabel                            J        e/f
LocalQNFormat                           %5.1f    %2s
                           

# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
None


# Calculate stick spectra or cross sections #
Temperature                             300
Range                                   0          30000
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.001          # If Y, default value 0.001
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   N


# Calculate cross sections #
Pressure                                1
Npoints/BinSize                         Npoints    10001
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt        
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
Cutoff(Y/N)                             Y          100            # If Y, default value 25
QNsFilter(Y/N)                          N          par[+]   e/f[e]   v[0,1,2,3]  
DopplerHWHM(Y/N)                        Y          0.1            # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   N
```