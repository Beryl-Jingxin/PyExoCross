# Examples

The examples of the whole input files for the ExoMol and HITRAN databases.

## Example for the ExoMol database

```bash
# Data source #
Database                                ExoMol
Molecule                                H2O
Isotopologue                            1H2-16O
Dataset                                 POKAZATEL
MolIsoID                                11


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           1
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
Non-LTE                                 0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               32
NCPUfiles                               32
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                Ka      Kc      v1      v2      v3      Gamma_rve
QNsformat                               %2d     %2d     %2d     %2d     %2d     %2s


# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    10000                     # Maximal temperature in K 


# Calculate lifetimes #
Compress(Y/N)                           Y                         # If Y, save as .states.bz2 file; otherwise, save as .states file


# Calculate oscillator strengths #
gf/f                                    gf
PlotOscillatorStrength(Y/N)             N  
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30


# Calculate stick spectra or cross sections #
Temperature                             300
Range                                   0          41200          # Unit cm-1
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          N          Ka[]  Kc[]  v1[]  v2[1,;,0]  v3[]  Gamma_rve[]


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30


# Calculate non-LTE stick spectra #
Tvib                                   2000
Trot                                   296
QNsVibLabel                            v,eS
QNsRotLabel                            J,e/f  
PlotNonLTE(Y/N)                        Y
Y-axisLimitNonLTE                      1e-30                     # Default value is 1e-30


# Calculate cross sections #
Pressure                                1
Npoints/BinSize                         BinSize   0.1
Broadeners                              H2       He  
Ratios                                  0.75     0.15  
Profile                                 SciPyVoigt   
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25, unit cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   N
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30

```

## Example for the HITRAN database

```bash
# Data source #
Database                                HITRAN
Molecule                                H2O
Isotopologue                            1H2-16O
Dataset                                 H2O-HITRAN
MolIsoID                                11


# File path #
ReadPath                                /home/jingxin/data/HITRAN/H2O.par
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
Non-LTE                                 0
CrossSections                           1


# Cores and chunks #
NCPUtrans                               32
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                v1     v2    v3     J     Ka    Kc    F    Sym  
QNsformat                               %2d    %2d   %2d    %3d   %3d   %3d   %5s  %1s


# Conversion #  
ConversionFormat                        2  
ConversionFrequncyRange                 0          1000   
GlobalQNLabel                           v1     v2    v3   
GlobalQNFormat                          %2d    %2d   %2d   
LocalQNLabel                            J     Ka    Kc    F    Sym  
LocalQNFormat                           %3d   %3d   %3d   %5s  %1s   
ConvUncFilter(Y/N)                      N          0.005          # If Y, default value 0.01
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30   
  

# Calculate partition functions, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
Compress(Y/N)                           N                         # If Y, save as .states.bz2 file; otherwise, save as .states file


# Calculate oscillator strengths #
gf/f                                    gf
PlotOscillatorStrength(Y/N)             Y  
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30


# Calculate stick spectra or cross sections #
Temperature                             300
Range                                   0          1000           # Unit cm-1
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          Y          v1[]       v2[]       v3[1,0;2,]


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y
Y-axisLimitStickSpectra                 1e-40                     # Default value is 1e-30


# Calculate non-LTE stick spectra #
Tvib                                   2000
Trot                                   296
QNsVibLabel                            v,eS
QNsRotLabel                            J,e/f  
PlotNonLTE(Y/N)                        Y
Y-axisLimitNonLTE                      1e-30                     # Default value is 1e-30


# Calculate cross sections #
Pressure                                1
Npoints/BinSize                         Npoints    10001
Broadeners                              Air        Self 
Ratios                                  0.7        0.3   
Profile                                 Lorentzian  
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
Cutoff(Y/N)                             Y          25             # If Y, default value 25, unit cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30
```
