# Non-LTE

`NLTEMethod` selects the population model:

`T`: two-temperature distribution using `Tvib` and `Trot`.

`D`: custom vibrational density $n_{\textrm{vib}}$ with a rotational
temperature `Trot`. Give the density-file path after `D`. `Temperatures` and
`Tvib` are not used.

`P`: custom rovibrational populations. Give the population-file path after
`P`. `Tvib` and `Trot` are not used.

``QNsVibLabel``: Please provide vibrational quantum number labels seperated by ``,``.

``QNsRotLabel``: Please provide rotational quantum number labels seperated by ``,``.

For ExoMol and ExoAtom input, these labels use the metadata-derived short names, for example `Herzberg:n1` should be written as `n1`. The `QNslabel` and `QNsformat` rows in examples are optional overrides for custom or legacy files.

## Two-temperature grids

`Tvib`: vibrational temperatures in K.

`Trot`: rotational temperatures in K.

| `Tvib` and `Trot`  |          | T = ?                                  |
| :----------------- | -------- | -------------------------------------- |
| 296                | &#x2705; | 296 K                                  |
| 300,1000,3000,8000 | &#x2705; | 300 K, 1000 K, 3000 K, 8000 K          |
| 1000:5000:1000     | &#x2705; | 1000 K, 2000 K, 3000 K, 4000 K, 5000 K |
| 300, 3000          | &#x274C; |                                        |
| 1000: 3000: 500    | &#x274C; |                                        |
| [300,3000]         | &#x274C; |                                        |
| [300, 3000]        | &#x274C; |                                        |

$N_{T_\textrm{vib}}$ and $N_{T_\textrm{rot}}$ and $N_P$ can be different. The total number of result files should be: 

$$N_\textrm{Total} = N_{T_\textrm{vib}} \times N_{T_\textrm{rot}} \times N_P$$ 

Here $N_P$ is the number of pressures if required.

## Two temperatures Treanor distribution

The state energy is the sum of the rotational and vibrational state energy:

$$
    \tilde{E}^{\textrm{tot}} = \tilde{E}^{\textrm{vib}}_{\textrm{QN}_{\textrm{vib}}} + \tilde{E}^{\textrm{rot}}_{\textrm{QN}_{\textrm{rot}}}.
$$


The non-LTE partition function equation is:

$$
    Q(T) = \sum_n g_n^{\textrm{tot}} e^{-c_2\tilde{E}^{\textrm{vib}}_{\textrm{QN}_{\textrm{vib}}}/T_{\textrm{vib}}} 
    e^{-c_2\tilde{E}^{\textrm{rot}}_{\textrm{QN}_{\textrm{rot}}}/T_{\textrm{rot}}}. 
$$


For absorption, define the vibrational and rotational transition energies as
$\Delta\tilde{E}_{\textrm{vib}}=\tilde{E}'_{\textrm{vib}}-\tilde{E}''_{\textrm{vib}}$
and
$\Delta\tilde{E}_{\textrm{rot}}=\tilde{E}'_{\textrm{rot}}-\tilde{E}''_{\textrm{rot}}$.
The intensity is:

$$
    I(f \gets i) = \frac{g'{A}_{fi}}{8 \pi c \tilde{v}^2_{fi}}
    \frac{e^{-c_2(\tilde{E}_{\textrm{vib}}''/T_{\textrm{vib}}
    +\tilde{E}_{\textrm{rot}}''/T_{\textrm{rot}})}}{Q(T)}
    \left[1-e^{-c_2(\Delta\tilde{E}_{\textrm{vib}}/T_{\textrm{vib}}
    +\Delta\tilde{E}_{\textrm{rot}}/T_{\textrm{rot}})}\right].
$$

The emissivity equation is:

$$
    \varepsilon (i \gets f) = \frac{g'{A}_{fi}hc}{4 \pi} 
    \frac{e^{-c_2 \tilde{E}_{\textrm{rot}}' / T_{\textrm{rot}}} e^{-c_2 \tilde{E}_{\textrm{vib}}' / T_{\textrm{vib}}}}{Q(T)}.
$$

## Pressure broadening for non-LTE cross sections

The population model and pressure-broadening model are independent. Methods
`T`, `D`, and `P` use the same broadening inputs as LTE cross sections. For
example, fixed reference parameters for an H$_2$/He mixture can be supplied as:

```bash
Broadeners                              H2       He
Ratios                                  0.85     0.15
GammaBroadening                         0.06     0.03
NBroadening                             0.60     0.30
Profile                                 SciPyVoigt
LorentzianHWHM(Y/N)                     N          0.5
```

The four lists correspond by position. `GammaBroadening` is the reference
Lorentzian HWHM coefficient in cm$^{-1}$ bar$^{-1}$ at 296 K;
`NBroadening` is dimensionless. For each spectrum condition, PyExoCross uses
the active profile temperature and pressure to evaluate

$$
\gamma(T,P)=\sum_i x_i\gamma_{0,i}
\left(\frac{296\ \mathrm{K}}{T}\right)^{n_i}
\frac{P}{1\ \mathrm{bar}}.
$$

For method `T`, the profile temperature is `Trot`. Method `D` also uses
`Trot`. Custom broadening is ignored by Doppler and Gaussian profiles. Setting
`LorentzianHWHM(Y/N)` to `Y` supplies a final constant HWHM and takes priority.

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
SpeciesID                               501


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/MgH_ExoMol_nlte.log    # Use None to disable log-file output


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
Device                                  CPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default), CUDA, PyTorch-CUDA, CuPy-CUDA, or MPS (used only when Device=GPU)
GPUBatchLines                           8192                      # GPU line-batch size (only used when Device=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when Device=GPU)


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                +/-  e/f   ElecState   v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %12s        %3d   %3d      %5.1f    %5.1f


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             Non-LTE                   # 'LTE' or 'Non-LTE'
Temperatures                            300,3000                  # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          e/f[]   v[0,;1,;2,;3,;4,;,0;,1;,2;,3;,4] 


# Calculate non-LTE #
NLTEMethod                              T                         # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                    300,3000,5000
Trot                                    300,500
QNsVibLabel                             v,ElecState
QNsRotLabel                             J,e/f   


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule


# Calculate cross sections #
Pressures                               1,10                      # Pressures in unit bar
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
Database                                ExoMolHR
Molecule                                NO
Isotopologue                            14N-16O
SpeciesID                               81


# File path #
ReadPath                                /mnt/data/exomolhr/exomolhr_results/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross//log/NO_ExoMolHR_nlte.log    # Use None to disable log-file output


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
Device                                  CPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default), CUDA, PyTorch-CUDA, CuPy-CUDA, or MPS (used only when Device=GPU)
GPUBatchLines                           8192                      # GPU line-batch size (only used when Device=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when Device=GPU)


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            296,1000                  # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   24         53452          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          


# Calculate non-LTE #
NLTEMethod                              T                         # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                    1000,2000
Trot                                    296
QNsVibLabel                             v,ElecState
QNsRotLabel                             J,e/f            


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule


# Calculate cross sections #
Pressures                               1                         # Pressure in unit bar
Npoints/BinSize                         BinSize    1              # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  0.1        
Profile                                 SciPyVoigt      
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             N          25             # If Y, default value 25 cm-1
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
SpeciesID                               81


# File path #
ReadPath                                /home/jingxin/data/HITRAN/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_HITRAN_nlte.log    # Use None to disable log-file output



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
Device                                  CPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default), CUDA, PyTorch-CUDA, CuPy-CUDA, or MPS (used only when Device=GPU)
GPUBatchLines                           8192                      # GPU line-batch size (only used when Device=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when Device=GPU)


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                J       X     Omega   v1      Sym    F
QNsformat                               %5.1f   %2s   %3s     %2d     %1s    %5s


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             Non-LTE                   # 'LTE' or 'Non-LTE'
Temperatures                            1000,2000                 # Temperature in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          N          X[]       Omega[]       v1[1,0;2,]


# Calculate non-LTE #
NLTEMethod                              T                         # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Tvib                                    1000:2000:1000
Trot                                    300,400
QNsVibLabel                             v1,X
QNsRotLabel                             J    


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule



# Calculate cross sections #
Pressures                               1                         # Pressures in unit bar
Npoints/BinSize                         BinSize    0.1            # Same unit as WnWlUnit
Broadeners                              Air        Self
Ratios                                  0.7        0.3
Profile                                 SciPyVoigt        
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        Y          3              # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     N          0.5            # Calculate pressure-dependent Lorentzian HWHM
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

For absorption, the rotational population and stimulated-emission correction
both use `Trot`:

$$
    I(f \gets i) = \frac{g'{A}_{fi}}{8\pi c\tilde{v}_{fi}^{2}}
    \frac{n_{\textrm{vib}}''e^{-c_2\tilde{E}_{\textrm{rot}}''/T_{\textrm{rot}}}}{Q}
    \left(1-e^{-c_2\tilde{v}_{fi}/T_{\textrm{rot}}}\right).
$$

Only `Trot` determines the temperature grid in method `D`.

*Example*

```bash
# Data source #
Database                                ExoAtom
Atom                                    Ar
Dataset                                 NIST


# File path #
ReadPath                                /home/jingxin/data/NLTE/         #/mnt/data/exoatom/exoatom_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/Ar_QDB_nlte_density.log    # Use None to disable log-file output


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
Device                                  CPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default), CUDA, PyTorch-CUDA, CuPy-CUDA, or MPS (used only when Device=GPU)
GPUBatchLines                           8192                      # GPU line-batch size (only used when Device=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when Device=GPU)


# Quantum numbers for conversion, stick spectra, and cross sections#
QNslabel                                configuration     Multiple     parity
QNsformat                               %50s              %30s         %2s  


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             Non-LTE                   # LTE or Non-LTE
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   8000       15500          # Same unit as WnWlUnit
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          Y          configuration[]     Multiple[]     parity[]


# Calculate non-LTE #
NLTEMethod                              D           /home/jingxin/data/NLTE/Ar/NIST/Ar_density.csv     # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
Trot                                    296
QNsVibLabel                             nan
QNsRotLabel                             nan       


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule


# Calculate cross sections #
Pressures                               1                         # Pressures in unit bar
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

Method `P` reads the supplied state populations directly. `Temperatures` is
retained for the requested spectrum conditions, while `Tvib` and `Trot` are
not used by the population model.

*Example*

```bash
# Data source #
Database                                ExoAtom
Atom                                    Ar
Dataset                                 NIST


# File path #
ReadPath                                /home/jingxin/data/NLTE/         #/mnt/data/exoatom/exoatom_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/Ar_QDB_nlte_pop.log    # Use None to disable log-file output


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
Device                                  CPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default), CUDA, PyTorch-CUDA, CuPy-CUDA, or MPS (used only when Device=GPU)
GPUBatchLines                           8192                      # GPU line-batch size (only used when Device=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when Device=GPU)


# Quantum numbers for conversion, stick spectra, and cross sections#
QNslabel                                configuration     Multiple     parity
QNsformat                               %50s              %30s         %2s  


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             Non-LTE                   # LTE or Non-LTE
Temperatures                            2000                      # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   8000       15500          # Same unit as WnWlUnit
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          Y          configuration[]     Multiple[]     parity[]


# Calculate non-LTE #
NLTEMethod                              P           /home/jingxin/data/NLTE/Ar/NIST/Ar_Ids.csv     # 'T'(TvibTrot) or 'D'(Density) or 'P'(Population)
QNsVibLabel                             nan
QNsRotLabel                             nan


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule


# Calculate cross sections #
Pressures                               1                         # Pressures in unit bar
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
