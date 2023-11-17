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
Conversion                              1
PartitionFunctions                      1
SpecificHeats                           1
CoolingFunctions                        1
Lifetimes                               1
OscillatorStrengths                     1
StickSpectra                            1
CrossSections                           1


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f


# Conversion #
ConversionFormat                        1  
ConversionFrequncyRange                 0          30000  
GlobalQNLabel                           eS       v        Omega
GlobalQNFormat                          %9s     %2d      %4s
LocalQNLabel                            J        e/f
LocalQNFormat                           %5.1f    %2s
ConvUncFilter(Y/N)                      Y          0.01           # If Y, default value 0.01
ConvThreshold(Y/N)                      Y          1e-30          # If Y, default value 1e-30
   

# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
Compress(Y/N)                           Y                         # If Y, save as .states.bz2 file; otherwise, save as .states file


# Calculate oscillator strengths #
fg/f                                    fg
Ncolumns                                4                         # 3 (without A) or 4 (with A)


# Calculate stick spectra or cross sections #
Temperature                             300
Range                                   0          30000
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          Y          par[]   e/f[e,f]   v[1,;2,;,0;4,4;4,3]  


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y
Y-axisLimitStick                        1e-30                     # Default value is 1e-30


# Calculate cross sections #
Pressure                                1
Npoints/BinSize                         Npoints    10001
Broadeners                              Default  
Ratios                                  1.0  
Profile                                 SciPyVoigt  
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
Cutoff(Y/N)                             Y          100            # If Y, default value 25
DopplerHWHM(Y/N)                        Y          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y
Y-axisLimitXsec                         1e-40                     # Default value is 1e-30
```

## Example for the HITRAN database

```bash
# Data source #
Database                                HITRAN
Molecule                                H2S
Isotopologue                            1H2-32S
Dataset                                 H2S-HITRAN
MolIsoID                                311


# File path #
ReadPath                                /home/jingxin/data/HITRAN/H2S.par
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            1
CrossSections                           1


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                v1     v2    v3     Ka     Kc    F     Sym
QNsformat                               %2d    %2d   %2d    %3d    %3d   %5s   %1s  


# Conversion #
ConversionFormat                        2  
ConversionFrequncyRange                 0         12000  
GlobalQNLabel                           v1       v2       v3
GlobalQNFormat                          %2d      %2d      %2d
LocalQNLabel                            J        Ka       Kc      F      Sym
LocalQNFormat                           %3d      %3d      %3d     %5s    %1s
ConvUncFilter(Y/N)                      N          0.005          # If Y, default value 0.01
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30  


# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    3000                      # Maximal temperature in K 


# Calculate lifetimes #
Compress(Y/N)                           Y                         # If Y, save as .states.bz2 file; otherwise, save as .states file


# Calculate oscillator strengths #
fg/f                                    fg
Ncolumns                                4                         # 3 (without A) or 4 (with A)


# Calculate stick spectra or cross sections #
Temperature                             300
Range                                   0          12000
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          Y          v1[1,]   v2[1,]   v3[1,0;1,1] 


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y
Y-axisLimitStick                        1e-40                     # Default value is 1e-30


# Calculate cross sections #
Pressure                                1
Npoints/BinSize                         Npoints    10001
Broadeners                              Default
Ratios                                  1.0
Profile                                 Gaussian  
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
Cutoff(Y/N)                             Y          25             # If Y, default value 25
DopplerHWHM(Y/N)                        Y          0.321          # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   N
Y-axisLimitXsec                                                   # Default value is 1e-30
```
