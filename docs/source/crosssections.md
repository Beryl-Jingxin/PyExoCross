Cross sections
==============

`Temperature`: Please provide temperature in unit K.

`Range`: Give two values as the minimum and maximum of the wavenumber range in unit cm-1. No `,` or `;` between these two numbers, just leave blank here.

`Absorption/Emission`: Choose `Absorption` or `Emission`.

`Pressure`: Please provide pressure in unit bar.

`Npoints/BinSize`: `Npoints` is the number of the points in grid. `BinSiza` is the interval size of the grid.

`Wavenumber(wn)/wavelength(wl)`: Choose `wn` or `wl`.

If  `PredissocXsec(Y/N)` is yes, predissociation lifetimes will be used/calculated when calculating cross sections with Voigt profile.

If  `Cutoff(Y/N)` is yes, you can provide cutoff here in unit cm-1.

If you want a figure of corss sections, please set `Y` for `PlotCrossSection(Y/N)`.

And if you want set the lower limit of y-axis for plotting, please write after `Y-axisLimitXsec`, otherwise, the default lower limit y-axis is 1e-30.

*Example*

```
# Calculate LTE or Non-LTE stick spectra or cross sections #
Temperature                             2000
Range                                   0          30000          # Unit cm-1
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          N          v[0,;1,;2,;3,;4,;,0;,1;,2;,3;,4] 


# Calculate cross sections #
Pressure                                1
Npoints/BinSize                         BinSize   0.1
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 SciPyVoigt        
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
PredissocXsec(Y/N)                      N
Cutoff(Y/N)                             Y          25             # If Y, default value 25, unit cm-1 
DopplerHWHM(Y/N)                        N          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y
Y-axisLimitXsec                         1e-40                     # Default value is 1e-30
```

## Filters

`Y/N`: `Y`, `YES`, `Yes`, `yes` and `N`, `NO`, `No`, `no` all can work. If you don't use it, write `N` here. You don't need to change the content behind it.
If after using filters, the program gets an empty result, then you will receive an warning to ask you write new filter values to icrease the range.

If `UncFilter(Y/N)` is yes, the value is the maximum uncertainty you require.

If `Threshold(Y/N)` is yes, the value is the minimum intensity you require.

If `QNsFilter(Y/N)` is yes, program will do filter on quantum numbers.
Write the quantum number labels required here, the spelling of the quantum number labels must be the same as `QNslabel`.
The other quantum number labels which are in `QNslabel` but not in `QNsFilter(Y/N)` will not be stored in the result file.
Write the quantum number values required after each label in `[]` and seperated by `,` and `;`, don't leave blank between different values inside the `[]`.
Leave blank between different quantum number labels, don't write `,`.
If you need all values of a quantum number label, write this label and wirte nothing inside the `[]`. Note, don't write any blank inside `[]`, you should write `[]`, not `[ ]`.
Inside `[]` use `,` and `;` don't write any blank inside the `[]`. Outside `[]`, use blank , don't write any `,` or `;` outside `[]`.
For one quantum number label, write in one `[]`, you can provide the quantum number values for upper and lower states, and seperated by `,`.
For one quantum number label, write in one `[]`, you can provide more than one pair of values, and seperated by `;`.

*Example*

`v1[]` means you want quantum number label v1 and you want all quantum number values of this label v1.`v1[1,0]` means you want quantum number label v1 and you want the upper QN = 1 and lower QN = 0. So v1' = 1 and v1" = 0.`v1[,0]` means you want quantum number label v1 and you want all upper QN but the lower QN = 0. So v1' = 0, 1, 2, 3, 4, ... and v1" = 0. `v1[3,]` means you want quantum number label v1 and you want all lower QN but the upper QN = 3. So v1' = 3 and v1" = 0, 1, 2, 3, 4, ... `v1[1,1;2,2]` means you want quantum number label v1 and you want when v1' = 1, v1" = 1; when v1' = 2, v1" = 2.`v1[1,;,0;5,5]  v2[]` means you want quantumnumber labels v1 and v2. For v1, you want all lines with v1' = 1 , all lines with v1" = 0 and the lines with v1' = 5 and at the same time v1" = 5. Meanwhile, you want all lines for v2.

* The definition file `.def` of ExoMol database (available at [exomol.com](https://www.exomol.com/)) provides the labels and formats of the quantum numbers for each species for reference.
* HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.

**Note**

You can define the quantum number column name by yourself, but please make sure it has letters without any spaces.
e.g. 'c1', 'c2', 'v1', 'v2', 'electronicState', 'electronic_state', '1v', '2v', 'M/E/C'.
Wrong format of the quantum number column nams: '1', '2', 'electronic state'.

*Example*

```bash
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          N          par[]   e/f[]   v[1,;2,2;2,1;,0]  
```

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

*Example*

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
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
Non-LTE                                 0
CrossSections                           1


# Cores and chunks #
NCPUtrans                               32
NCPUfiles                               32
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                Ka      Kc      v1      v2      v3      Gamma_rve
QNsformat                               %2d     %2d     %2d     %2d     %2d     %2s


# Calculate LTE or Non-LTE stick spectra or cross sections #
Temperature                             300
Range                                   0          41200          # Unit cm-1
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          N          Ka[]  Kc[]  v1[]  v2[1,;,0]  v3[]  Gamma_rve[]


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
QNslabel                                J       X     Omega   v1      m      Sym  
QNsformat                               %5s     %2s   %3s     %2d     %1s 


# Calculate LTE or Non-LTE stick spectra or cross sections #
Temperature                             1000
Range                                   1000       5000          # Unit cm-1
Absorption/Emission                     Emission                 # 'Absorption' or 'Emission'
UncFilter(Y/N)                          No         0.01          # If Y, default value 0.01
Threshold(Y/N)                          NO         1e-30         # If Y, default value 1e-30
QNsFilter(Y/N)                          N          par[]   e/f[e,e]   v[1,;2,;,0;4,4;4,3]  


# Calculate cross sections #
Pressure                                0.1
Npoints/BinSize                         BinSize    0.1
Broadeners                              Air        Self  
Ratios                                  0.7        0.3   
Profile                                 SciPyVoigt  
Wavenumber(wn)/wavelength(wl)           wl                        # 'wn' or 'wl'
PredissocXsec(Y/N)                      no
Cutoff(Y/N)                             N          100            # If Y, default value 25, unit cm-1
DopplerHWHM(Y/N)                        n          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     n          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   NO
Y-axisLimitXsec                                                   # Default value is 1e-30
```
