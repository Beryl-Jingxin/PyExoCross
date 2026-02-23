# Stick spectra

`LTE/Non-LTE`: Choose `LTE` or `Non-LTE`. \
If you choose `LTE`, please ignore `# Calculate non-LTE #` section. \
If you choose `Non-LTE`, more details can be found from [**Non-LTE**](`https://pyexocross.readthedocs.io/en/latest/nonlte.html`).

`WnWlUnit`: Choose to provide the range of wavenumber `wn` in unit of `cm-1` (cm⁻¹), or wavelength in unit of `um` (μm) or `nm`.

`Range`: Give two values as the minimum and maximum of the wavenumber range in unit of cm⁻¹ or wavelength range in unit of μm or nm. Please use the same unit as `WnWlUnit`. Don't use `,` or `;` between these two numbers, just leave blank here. 

`Absorption/Emission`: Choose `Absorption` or `Emission`.

The LTE intensity equation is:

$$
    I(f \gets i) = \frac{g'{A}_{fi}}{8 \pi c \tilde{v}^2_{fi}} 
    \frac{e^{-c_2 \tilde{E}'' / T} (1 - e^{-c_2 \tilde{v}_{fi} 
    / T })}{Q(T)}.
$$

The LTE emissivity equation is:

$$
    \varepsilon (i \gets f) = \frac{g'{A}_{fi}hc}{4 \pi}\frac{e^{-c_2 \tilde{E}'/T}}{Q(T)}.
$$

## Temperatures

`Temperatures`: Please provide temperatures in unit K.

| `Temperatures`     |          | T = ?                                  |
| :----------------- | -------- | -------------------------------------- |
| 296                | &#x2705; | 296 K                                  |
| 300,1000,3000,8000 | &#x2705; | 300 K, 1000 K, 3000 K, 8000 K          |
| 1000:5000:1000     | &#x2705; | 1000 K, 2000 K, 3000 K, 4000 K, 5000 K |
| 300, 3000          | &#x274C; |                                        |
| 1000: 3000: 500    | &#x274C; |                                        |
| [300,3000]         | &#x274C; |                                        |
| [300, 3000]        | &#x274C; |                                        |

## Filters

`UncFilter(Y/N)`, `Threshold(Y/N)`, and `QNsFilter(Y/N)` are filters, please see [**Filters**](`https://pyexocross.readthedocs.io/en/latest/filters.html`).

## Plot stick spectra

`PlotStickSpectra(Y/N)`: If you need plot a figure of stick spectra (intensity or emisivity), please write `Y` after `PlotStickSpectra(Y/N)`, otherwise, write `N`.

`PlotStickSpectraMethod`: Choose to plot y-axis in linear `lin` or logarithm `log`.

`PlotStickSpectraWnWl`: Choose to plot x-axis with wavenumber `wn` in unit of `cm-1` (cm⁻¹), or wavelength `wl` in unit of `um` (μm) or `nm`.

`Y-axisLimitStickSpectra`: If you need set the lower limit of y-axis for plotting, please write after `Y-axisLimitStickSpectra`, otherwise, the default lower limit y-axis is 1e-30 ($=10^{-30}$) in unit of cm/molecule.

*Example*

Save wavenumber in unit of cm⁻¹ in the file and 
plot logarithm stick spectra in unit of cm/molecule and wavenumber in unit of cm⁻¹.

```bash
# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            300,3000                  # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          par[]   e/f[]  eS[]  v[1,1;1,0;2,;,0]  


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule
```

Save wavenumber in unit of cm⁻¹ in the file and 
plot logarithm stick spectra in unit of cm/molecule and wavelength in unit of nm.

```bash
# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            1000:4000:1000            # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          Y          par[]   e/f[]  eS[]  v[]  


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wl         nm             # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule
```

Save wavenumber in unit of μm in the file and 
plot linear stick spectra in unit of cm/molecule and wavelength in unit of nm.

```bash
# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            300                       # Temperatures in unit of K
WnWlUnit                                wl         um             # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          Y          par[]   e/f[]  eS[]  v[]  


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  lin                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wl         nm             # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule
```

*Example*

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
LogFilePath                             /home/jingxin/data/pyexocross/log/MgH_ExoMol_stick_Ab.log


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            1
CrossSections                           0


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            300,3000                  # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   0          30000          # Same unit as WnWlUnit
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          par[]   e/f[]  eS[]  v[1,1;1,0;2,;,0]  


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule
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
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_HITRAN_stick_Em_T2000_wl150-1000.log


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            1
CrossSections                           0


# Cores and chunks #
NCPUtrans                               8
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                +/-    e/f     State    v     Lambda   Sigma    Omega
QNsformat                               %1s    %1s     %5s      %5d   %5d      %5.1f    %5.1f


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            2000:4000:500             # Temperatures in unit of K
WnWlUnit                                wl         nm             # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   150        1000           # Same unit as WnWlUnit
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01 cm-1
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          Y          State[]  v[]  Lambda[]  Sigma[]  Omega[]


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule
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
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_HITRAN_stick_Em_T1000_wn1000-5000.log


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            1
CrossSections                           0


# Cores and chunks #
NCPUtrans                               32
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                J       X     Omega   v1      m      Sym    
QNsformat                               %5s     %2s   %3s     %2d     %1s    %1s


# Calculate stick spectra or cross sections #
LTE/Non-LTE                             LTE                       # 'LTE' or 'Non-LTE'
Temperatures                            1000                      # Temperatures in unit of K
WnWlUnit                                wn         cm-1           # Wavenumber (wn in unit of cm-1) or wavelength (wl in unit of um or nm)
Range                                   1000       5000           # Same unit as WnWlUnit
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.005          # If Y, default value 0.01 cm-1
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30 cm/molecule
QNsFilter(Y/N)                          N          par[]   e/f[e,e]   v[1,1;1,0]  


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   N          
PlotStickSpectraMethod                  log                       # Plot in linear (lin) or logarithm (log)
PlotStickSpectraWnWl                    wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[nm or um])
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30 cm/molecule
```