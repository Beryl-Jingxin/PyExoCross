# Cross sections

`LTE/Non-LTE`: Choose `LTE` or `Non-LTE`. \
If you choose `LTE`, please ignore `# Calculate non-LTE #` section. \
If you choose `Non-LTE`, more details can be found from [**Non-LTE**](`https://pyexocross.readthedocs.io/en/latest/nonlte.html`).

`WnWlUnit`: Choose to provide the range of wavenumber `wn` in unit of `cm-1` (cm⁻¹), or wavelength in unit of `um` (μm) or `nm`.

`Range`: Give two values as the minimum and maximum of the wavenumber range in unit of cm⁻¹ or wavelength range in unit of μm or nm. Please use the same unit as `WnWlUnit`. Don't use `,` or `;` between these two numbers, just leave blank here.

`Absorption/Emission`: Choose `Absorption` or `Emission`.



`Npoints/BinSize`: `Npoints` is the number of the points in grid. `BinSiza` is the interval size of the grid , use the same unit as `WnWlUnit`.

`PredissocXsec(Y/N)`: If `PredissocXsec(Y/N)` is yes, predissociation lifetimes will be used or calculated when calculating cross sections with Voigt profile.

`Cutoff(Y/N)`: If `Cutoff(Y/N)` is yes, you can provide wing cutoff here in unit of cm⁻¹.

## Temperatures and pressures

`Temperatures`: Please provide temperatures in unit of K. 

| `Temperatures`     |          | T = ?                                  |
| :----------------- | -------- | -------------------------------------- |
| 296                | &#x2705; | 296 K                                  |
| 300,1000,3000,8000 | &#x2705; | 300 K, 1000 K, 3000 K, 8000 K          |
| 1000:5000:1000     | &#x2705; | 1000 K, 2000 K, 3000 K, 4000 K, 5000 K |
| 300, 3000          | &#x274C; |                                        |
| 1000: 3000: 500    | &#x274C; |                                        |
| [300,3000]         | &#x274C; |                                        |
| [300, 3000]        | &#x274C; |                                        |

`Pressures`: Please provide pressures in unit bar.

| `Pressures`      |          | P = ?                               |
| :--------------- | -------- | ----------------------------------- |
| 1                | &#x2705; | 1 bar                               |
| 1.0              | &#x2705; | 1.0 bar                             |
| 0.5,1.0,5.0,10.0 | &#x2705; | 0.5 bar, 1.0 bar, 5.0 bar, 10.0 bar |
| 1:5:1            | &#x2705; | 1 bar, 2 bar, 3 bar, 4 bar, 5 bar   |
| 1, 10            | &#x274C; |                                     |
| 1: 5: 1          | &#x274C; |                                     |
| [1,10]           | &#x274C; |                                     |
| [1, 10]          | &#x274C; |                                     |

$N_T$ and $N_P$ can be different. The total number of result files should be:

$$N_\textrm{Total} = N_T \times N_P$$

## Filters

``UncFilter(Y/N)``, ``Threshold(Y/N)``, and ``QNsFilter(Y/N)`` are filters, please see [**Filters**](`https://pyexocross.readthedocs.io/en/latest/filters.html`).

## Broadeners

`Default` means:

1. Default value of temperature exponent for all lines is 0.5.
2. Default value of Lorentzian half-width for all lines is 0.07.

If you want to use the default values of the broadeners, please set the `Ratios` of `Broadeners` under `Default` as `1.0`.

The broadening types and ratio values are corresponding, please write them in order.

*Example*

```bash
Broadeners                              Default   
Ratios                                  1.0  
```

```bash
Broadeners                              Air      Self  
Ratios                                  0.7      0.3  
```

```bash
Broadeners                              H2       He   
Ratios                                  0.9      0.1   
```

## Line profiles

Choose line profile from:

`Doppler`, `Gaussian`, `Lorentzian`, `SciPyVoigt`, `SciPyWofzVoigt`, `HumlicekVoigt`, `ThompsonPseudoVoigt`, `KielkopfPseudoVoigt`, `OliveroPseudoVoigt`, `LiuLinPseudoVoigt`, `RoccoPseudoVoigt`, `BinnedDoppler`, `BinnedGaussian`, `BinnedLorentzian`, `BinnedVoigt`. Please note: no blank when you write the line profile name.

`DopplerHWHM(Y/N)`:

1. Doppler profile uses Doppler HWHM calculated by program, if you want to use Doppler profile, set `N` after `DopplerHWHM(Y/N)`. 
2. If you use Gaussian profile, please set `Y` or `U` after `DopplerHWHM(Y/N)`, your Doppler HWHM value will be used for calculating Gaussian profile.
3. If you prefer to provide a Doppler HWHM as a constant, write `Y` here and give the HWHM value. 
4. If you want to use your own Doppler HWHM list, please add them to the transitions file(s) as a new column on the right. Write `U` here and give this new column index (count start from 0). 
5. If you need the program to calculate tempertaure-dependent Doppler HWHM, please write `N` here.

`LorentzianHWHM(Y/N)`: 

1. If you prefer to provide a Lorentzian HWHM as a constant, write `Y` here and give the HWHM value. 
2. If you want to use your own Lorentzian HWHM list, please add them to the transitions file(s) as a new column on the right. Write `U` here and give this new column index (count start from 0). 
3. If you need the program to calculate tempertaure and pressure-dependent Loretzian HWHM, please write `N` here.

*Example*

```bash
Profile                                 SciPyVoigt
DopplerHWHM(Y/N)                        Y          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     Y          0.5            # Set Lorentzian HWHM as a constant
```

```bash
Profile                                 Doppler
DopplerHWHM(Y/N)                        N          0.3            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
```

```bash
Profile                                 Gaussian
DopplerHWHM(Y/N)                        U          4              # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
```

```bash
Profile                                 SciPyVoigt
DopplerHWHM(Y/N)                        U          4              # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     U          5              # Set Lorentzian HWHM as a constant
```

## Plot cross sections

`PlotCrossSection(Y/N)`: If you want a figure of corss sections, please set `Y`, otherwise, please write `N`.

`PlotCrossSectionMethod`: Choose to plot y-axis in linear `lin` or logarithm `log`.

`PlotCrossSectionWnWl`: Choose to plot x-axis with wavenumber `wn` in unit of `cm-1` (cm⁻¹), or wavelength `wl` in unit of `um` (μm) or `nm`.

`Y-axisLimitXsec`: If you want set the lower limit of y-axis for plotting, please write after `Y-axisLimitXsec`, otherwise, the default lower limit y-axis is 1e-30 ($=10^{-30}$) cm²/molecule.

*Example*

Save wavenumber in unit of cm⁻¹ in the file and 
plot logarithm cross sections in unit of cm²/molecule and wavenumber in unit of cm⁻¹.

```bash
# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            2000,5000                 # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          v[0,;1,;2,;3,;4,;,0;,1;,2;,3;,4] 


# Calculate cross sections #
Pressures                               0.1,1,10                  # Pressures in unit bar
Npoints/BinSize                         BinSize    0.1            # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt      
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-25                     # Default value is 1e-30 cm2/molecule
```

Save wavenumber in unit of cm⁻¹ in the file and 
plot logarithm cross sections in unit of cm²/molecule and wavelength in unit of nm.

```bash
# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            1000:3000:1000            # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          v[0,;1,;2,;3,;4,;,0;,1;,2;,3;,4] 


# Calculate cross sections #
Pressures                               1                         # Pressures in unit bar
Npoints/BinSize                         Npoint     30001          # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt      
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wl         nm             # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-25                     # Default value is 1e-30 cm2/molecule
```

Save wavelength in unit of μm in the file and 
plot linear cross sections in unit of cm²/molecule and wavelength in unit of nm.

```bash
# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            2000                      # Temperatures in unit of K
WnWlUnit                                wl         um             # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   6          10             # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          


# Calculate cross sections #
Pressures                               1                         # Pressure in unit bar
Npoints/BinSize                         BinSize    0.001          # Same unit as WnWlUnit
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt      
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  lin                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wl         nm             # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-25                     # Default value is 1e-30 cm2/molecule
```

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                H2O
Isotopologue                            1H2-16O
Dataset                                 POKAZATEL
SpeciesID                               11


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/H2O_ExoMol_xsec.log


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
NCPUfiles                               32
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                Ka      Kc      v1      v2      v3      Gamma_rve
QNsformat                               %2d     %2d     %2d     %2d     %2d     %2s


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            300                       # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          41200          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          Ka[]  Kc[]  v1[]  v2[1,;,0]  v3[]  Gamma_rve[]


# Calculate cross sections #
Pressures                               0.1,1                     # Pressures in unit bar
Npoints/BinSize                         BinSize    0.1            # Same unit as WnWlUnit
Broadeners                              H2       He  
Ratios                                  0.75     0.15  
Profile                                 SciPyVoigt   
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   N
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30 cm2/molecule
```

```bash
# Data source #
Database                                ExoMol
Molecule                                NO
Isotopologue                            14N-16O
Dataset                                 XABC
SpeciesID                               81


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_HITRAN_xsec_T2000_wl150-1000.log


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
NCPUtrans                               8
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                +/-    e/f     State    v     Lambda   Sigma    Omega
QNsformat                               %1s    %1s     %5s      %5d   %5d      %5.1f    %5.1f


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            300,2000,5000             # Temperatures in unit of K
WnWlUnit                                wl         nm             # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   150        1000           # Same unit as WnWlUnit
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          State[]  v[]  Lambda[]  Sigma[]  Omega[]


# Calculate cross sections #
Pressures                               1,10                      # Pressures in unit bar
Npoints/BinSize                         BinSize    0.1            # Same unit as WnWlUnit
Broadeners                              air    
Ratios                                  1.0        
Profile                                 Gaussian       
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          100            # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        Y          3              # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y          
PlotCrossSectionMethod                  lin                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wl         nm             # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
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
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_HITRAN_xsec_T1000_wn1000-5000.log


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
QNslabel                                J       X     Omega   v1      m      Sym  
QNsformat                               %5s     %2s   %3s     %2d     %1s 


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                      # 'LTE' or 'Non-LTE'
Temperatures                            1000:4000:1000           # Temperatures in unit of K
WnWlUnit                                wn         cm-1          # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   1000       5000          # Same unit as WnWlUnit
Absorption/Emission                     Emission                 # 'Absorption' or 'Emission'
UncFilter(Y/N)                          No         0.01          # If Y, default value 0.01 cm-1
Threshold(Y/N)                          NO         1e-30         # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          par[]   e/f[e,e]   v[1,;2,;,0;4,4;4,3]  


# Calculate cross sections #
Pressures                               0.1                       # Pressures in unit bar
Npoints/BinSize                         BinSize    0.1            # Same unit as WnWlUnit
Broadeners                              Air        Self  
Ratios                                  0.7        0.3   
Profile                                 SciPyVoigt  
PredissocXsec(Y/N)                      no
Cutoff(Y/N)                             N          100            # If Y, default value 25 cm-1
DopplerHWHM(Y/N)                        n          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     n          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   NO
PlotCrossSectionMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotCrossSectionWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitXsec                         1e-30                     # Default value is 1e-30 cm2/molecule
```
