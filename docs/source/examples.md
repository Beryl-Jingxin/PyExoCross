# Examples

The examples of the whole input files for databases can be obtained from the `input` folder of *PyExoCross*'s ***GitHub*** : [https://github.com/ExoMol/PyExoCross](https://github.com/ExoMol/PyExoCross).

## Example for the ExoMol database

```bash
# Data source #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
MolIsoID                                666  


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


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f


# Conversion #
ConversionFormat                        1  
ConversionFrequncyRange                 0          30000          # Wavenumber in unit of cm-1      
GlobalQNLabel                           eS      v      Omega
GlobalQNFormat                          %9s     %2d    %4s
LocalQNLabel                            J       e/f
LocalQNFormat                           %5.1f   %2s
ConvUncFilter(Y/N)                      N          0.01           # If Y, default value 0.01 cm-1
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30 cm/molecule
                           

# Calculate partition functions, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
Compress(Y/N)                           N                         # If Y, save as .states.bz2 file; otherwise, save as .states file


# Calculate oscillator strengths #
gf/f                                    f
PlotOscillatorStrength(Y/N)             Y         
PlotOscillatorStrengthMethod            log                       # Plot in linear (lin) or logarithm (log)
PlotOscillatorStrengthWnWl              wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[um or nm])
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperature                             2000                      # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          par[]  e/f  eS[]  v[0,;1,;2,;3,;4,;,0;,1;,2;,3;,4] 


# Calculate non-LTE #
NLTEMethod                             T                          # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                   2000
Trot                                   296
QNsVibLabel                            v,eS
QNsRotLabel                            J,e/f            


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule


# Calculate cross sections #
Pressure                                1                         # Pressure in unit bar
Npoints/BinSize                         BinSize    0.1            # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt      
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        Y          3              # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30 cm2/molecule
```

## Example for the ExoAtom database

```bash
# Data source #
Database                                ExoAtom
Atom                                    Li
Dataset                                 NIST


# File path #
ReadPath                                /mnt/data/exoatom/exoatom_data/
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


# Cores and chunks #
NCPUtrans                               1
NCPUfiles                               1
ChunkSize                               10000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                configuration     LS       parity  
QNsformat                               %30s              %30s     %2s      


# Conversion #
ConversionFormat                        1  
ConversionFrequncyRange                 0          43000          # Wavenumber in unit of cm-1      
GlobalQNLabel                           configuration     LS       
GlobalQNFormat                          %30s              %30s     
LocalQNLabel                            J       parity
LocalQNFormat                           %5.1f   %2s
ConvUncFilter(Y/N)                      N          0.01           # If Y, default value 0.01 cm-1
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30 cm/molecule


# Calculate partition functions, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
Compress(Y/N)                           N                         # If Y, save as .states.bz2 file; otherwise, save as .states file


# Calculate oscillator strengths #
gf/f                                    gf
PlotOscillatorStrength(Y/N)             Y         
PlotOscillatorStrengthMethod            log                       # Plot in linear (lin) or logarithm (log)
PlotOscillatorStrengthWnWl              wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[um or nm])
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperature                             2000                      # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          43000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N                   


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule


# Calculate cross sections #
Pressure                                1                         # Pressure in unit bar
Npoints/BinSize                         BinSize    0.1            # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 Gaussian      
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        Y          3              # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30 cm2/molecule
```


## Example for the HITRAN and HITEMP databases

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
ConversionFrequncyRange                 0          1000           # Wavenumber in unit of cm-1      
GlobalQNLabel                           v1     v2     v3   
GlobalQNFormat                          %2d    %2d    %2d   
LocalQNLabel                            J      Ka     Kc    F    Sym  
LocalQNFormat                           %3d    %3d    %3d   %5s  %1s   
ConvUncFilter(Y/N)                      N          0.005          # If Y, default value 0.01 cm-1
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30 cm/molecule
  

# Calculate partition functions, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
Compress(Y/N)                           N                         # If Y, save as .states.bz2 file; otherwise, save as .states file


# Calculate oscillator strengths #
gf/f                                    f
PlotOscillatorStrength(Y/N)             Y         
PlotOscillatorStrengthMethod            log                       # Plot in linear (lin) or logarithm (log)
PlotOscillatorStrengthWnWl              wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[um or nm])
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperature                             300                       # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          1000           # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          Y          v1[]    v2[]    v3[1,0;2,]             


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-35                     # Default value is 1e-30 cm/molecule


# Calculate non-LTE stick spectra #
NLTEMethod                             T                          # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                   2000
Trot                                   296
QNsVibLabel                            v1,v2,v3
QNsRotLabel                            J,Ka,Kc


# Calculate cross sections #
Pressure                                1                         # Pressure in unit bar
Npoints/BinSize                         Npoints    10001          # Same unit as WnWlUnit
Broadeners                              Air        Self 
Ratios                                  0.7        0.3         
Profile                                 Lorentzian      
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          3              # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30 cm2/molecule
```
