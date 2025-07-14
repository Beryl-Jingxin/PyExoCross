# Stick spectra

`LTE/Non-LTE`: Choose `LTE` or `Non-LTE`. \
If you choose `LTE`, please ignore `# Calculate non-LTE #` section. \
If you choose `Non-LTE`, more details can be found from [**Non-LTE**](`https://pyexocross.readthedocs.io/en/latest/nonlte.html`).

`Temperature`: Please provide temperature in unit K.

`Range`: Give two values as the minimum and maximum of the wavenumber range inunit $\textrm{cm}^{-1}$. Don't use `,` or `;` between these two numbers, just leave blank here.

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

## Filters

`UncFilter(Y/N)`, `Threshold(Y/N)`, and `QNsFilter(Y/N)` are filters, please see [**Filters**](`https://pyexocross.readthedocs.io/en/latest/filters.html`).

## Plot stick spectra

`PlotStickSpectra(Y/N)`: If you need plot a figure of stick spectra (intensity or emisivity), please write `Y` after `PlotStickSpectra(Y/N)`, otherwise, write `N`.

`PlotStickSpectra(Y/N)`: If you need set the lower limit of y-axis for plotting, please write after `Y-axisLimitStickSpectra`, otherwise, the default lower limit y-axis is 1e-30 ($=10^{-30}$) in unit cm/molecule.

*Example*

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
Temperature                             300
Range                                   0          30000
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.01           # If Y, default value 0.01
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          Y          par[]   e/f[]  eS[]  v[1,1;1,0;2,;,0]  


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   Y
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30
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
Temperature                             1000
Range                                   1000       5000
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          N          0.005          # If Y, default value 0.01
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          N          par[]   e/f[e,e]   v[1,1;1,0]  


# Calculate stick spectra #
PlotStickSpectra(Y/N)                   N
Y-axisLimitStickSpectra                 1e-30                     # Default value is 1e-30
```