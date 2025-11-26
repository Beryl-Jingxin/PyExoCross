# Non-LTE

`NLTEMethod`: 

`T`: Treanor distribution using two different temperatures `Tvib` and `Trot`. 

`D`: Vibronic distribution using `Trot` and custom vibrational density $n_{\textrm{vib}}$. After `D`, please give the custom vibrational density file path. Ignore `Tvib`.

`P`: Using custom rovibrational population. After `P`, please give the custom rovibrational population file path.  Ignore `Tvib` and `Trot`.

``Tvib``: Please provide vibrational temperature in unit of K.

``Trot``: Please provide rotational temperature in unit of K.

``QNsVibLabel``: Please provide vibrational quantum number labels seperated by ``,``.

``QNsRotLabel``: Please provide rotational quantum number labels seperated by ``,``.

## Two temperature Treanor distribution

The state energy is the sum of the rotational and vibrational state energy:

$$
    \tilde{E}^{\textrm{tot}} = \tilde{E}^{\textrm{vib}}_{\textrm{QN}_{\textrm{vib}}} + \tilde{E}^{\textrm{rot}}_{\textrm{QN}_{\textrm{rot}}}.
$$


The non-LTE partition function equation is:

$$
    Q(T) = \sum_n g_n^{\textrm{tot}} e^{-c_2\tilde{E}^{\textrm{vib}}_{\textrm{QN}_{\textrm{vib}}}/T_{\textrm{vib}}} 
    e^{-c_2\tilde{E}^{\textrm{rot}}_{\textrm{QN}_{\textrm{rot}}}/T_{\textrm{rot}}}. 
$$


The intensity equation is:

$$
    I(f \gets i) = \frac{g'{A}_{fi}}{8 \pi c \tilde{v}^2_{fi}} 
    \frac{e^{-c_2 \tilde{E}_{\textrm{rot}}'' / T_{\textrm{rot}}} e^{-c_2 \tilde{E}_{\textrm{vib}}'' / T_{\textrm{vib}}} (1 - e^{-c_2 \tilde{v}_{fi} / T_{\textrm{vib}} })}{Q(T)}.
$$

The emissivity equation is:

$$
    \varepsilon (i \gets f) = \frac{g'{A}_{fi}hc}{4 \pi} 
    \frac{e^{-c_2 \tilde{E}_{\textrm{rot}}' / T_{\textrm{rot}}} e^{-c_2 \tilde{E}_{\textrm{vib}}' / T_{\textrm{vib}}}}{Q(T)}.
$$

*Example*

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
LogFilePath                             /home/jingxin/data/pyexocross/log/MgH_ExoMol_nlte.log


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
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             Non-LTE                   # 'LTE' or 'Non-LTE'
Temperature                             300                       # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          e/f[]   v[0,;1,;2,;3,;4,;,0;,1;,2;,3;,4] 


# Calculate non-LTE #
NLTEMethod                              T                         # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                    2000
Trot                                    296
QNsVibLabel                             v,eS
QNsRotLabel                             J,e/f   


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

```bash
# Data source #
Database                                HITRAN
Molecule                                NO
Isotopologue                            14N-16O
Dataset                                 NO-HITRAN
MolIsoID                                81


# File path #
ReadPath                                /home/jingxin/data/HITRAN/NO.par
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_HITRAN_nlte.log



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
NCPUtrans                               32
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                J       X     Omega   v1      m      Sym    
QNsformat                               %5s     %2s   %3s     %2d     %1s    %1s


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             Non-LTE                   # 'LTE' or 'Non-LTE'
Temperature                             3000                      # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   1000       5000           # Same unit as WnWlUnit
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.005          # If Y, default value 0.01 cm-1
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          par[]   e/f[e,e]   v[1,1;1,0]  


# Calculate non-LTE #
NLTEMethod                              T                         # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                    2000
Trot                                    296
QNsVibLabel                             v,eS
QNsRotLabel                             J,e/f      


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule



# Calculate cross sections #
Pressure                                1                         # Pressure in unit bar
Npoints/BinSize                         Npoints    10001          # Same unit as WnWlUnit
Broadeners                              Air        Self 
Ratios                                  0.7        0.3     
Profile                                 Lorentzian        
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30 cm2/molecule
```

## Custom Vibrational density

The sum of the custom vibrational density $n_{\textrm{vib}}$ is normalized to one.

$$
    \sum_i n_i^{\textrm{vib}} = 1
$$

*Example*

```bash
# Data source #
Database                                ExoAtom
Atom                                    Ar
Dataset                                 NIST


# File path #
ReadPath                                /home/jingxin/data/NLTE/         #/mnt/data/exoatom/exoatom_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/Ar_QDB_nlte_density.log


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
NCPUtrans                               1
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra, and cross sections#
QNslabel                                configuration     Multiple     parity
QNsformat                               %50s              %30s         %2s  


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             Non-LTE                   # LTE or Non-LTE
Temperature                             2000                      # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   8000       15500          # Same unit as WnWlUnit
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          Y          configuration[]     Multiple[]     parity[]


# Calculate non-LTE #
NLTEMethod                              D           /home/jingxin/data/NLTE/Ar/NIST/Ar_density.csv     # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                    2000
Trot                                    296
QNsVibLabel                             v,eS
QNsRotLabel                             J,parity          


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule


# Calculate cross sections #
Pressure                                1                         # Pressure in unit bar
Npoints/BinSize                         BinSize    0.01           # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt    
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             N          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30 cm2/molecule
```

## Custom rovibrational population

*Example*

```bash
# Data source #
Database                                ExoAtom
Atom                                    Ar
Dataset                                 NIST


# File path #
ReadPath                                /home/jingxin/data/NLTE/         #/mnt/data/exoatom/exoatom_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/Ar_QDB_nlte_pop.log


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
NCPUtrans                               1
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra, and cross sections#
QNslabel                                configuration     Multiple     parity
QNsformat                               %50s              %30s         %2s  


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             Non-LTE                   # LTE or Non-LTE
Temperature                             2000                      # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   8000       15500          # Same unit as WnWlUnit
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          Y          configuration[]     Multiple[]     parity[]


# Calculate non-LTE #
NLTEMethod                              P           /home/jingxin/data/NLTE/Ar/NIST/Ar_Ids.csv     # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                    2000
Trot                                    296
QNsVibLabel                             v,eS
QNsRotLabel                             J,parity          


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule


# Calculate cross sections #
Pressure                                1                         # Pressure in unit bar
Npoints/BinSize                         BinSize    0.01           # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt       
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             N          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-40                     # Default value is 1e-30 cm2/molecule
```
