# Oscillator strengths

Weighted oscillator strength ``gf``, which is what people usually use.  
``f`` is the actual oscillator strength and the actual value of ``f`` is ``gf`` divided by ``g"``.
Sometimes people need oscillator strength ``f``, not ``gf``. 

The ressult file have 4 columns: 
ExoMol data: upper id, lower id, oscillator strength and wavenumber v.
HITRAN data: g', g", oscillator strength and wavenuber v.

``gf/f``: Choose weigthed oscillator strength ``gf`` or actual oscillator strength ``f``. 

## Plot oscillator strengths

``PlotOscillatorStrength(Y/N)``: If you need a oscillator strength figure, please write ``Y`` here. 

``PlotOscillatorStrengthMethod``: Choose to plot y-axis in linear `lin` or logarithm `log`.

``PlotOscillatorStrengthWnWl``: Choose to plot x-axis with wavenumber `wn` in unit of `cm-1` (cm⁻¹), or wavelength `wl` in unit of `um` (μm) or `nm`.

``Y-axisLimitOscillatorStrength``: If you want set the lower limit of y-axis for plotting, please write here, otherwise, the default lower limit y-axis is 1e-30.

*Example*

```bash
# Calculate oscillator strengths #
gf/f                                    gf
PlotOscillatorStrength(Y/N)             Y         
PlotOscillatorStrengthMethod            log                       # Plot in linear (lin) or logarithm (log)
PlotOscillatorStrengthWnWl              wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[um or nm])
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30
```

```bash
# Calculate oscillator strengths #
gf/f                                    f
PlotOscillatorStrength(Y/N)             Y         
PlotOscillatorStrengthMethod            lin                       # Plot in linear (lin) or logarithm (log)
PlotOscillatorStrengthWnWl              wl         um             # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[um or nm])
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30
```


The oscillator strengths equation is:

$$
    gf=\frac{g_\textrm{tot}' A_{fi}}{(c \tilde{v}_{fi})^2}.
$$

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                CO2
Isotopologue                            12C-16O2
Dataset                                 UCL-4000
MolIsoID                                21


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     1
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               4
ChunkSize                               1000000


# Calculate oscillator strengths #
gf/f                                    gf
PlotOscillatorStrength(Y/N)             Y         
PlotOscillatorStrengthMethod            log                       # Plot in linear (lin) or logarithm (log)
PlotOscillatorStrengthWnWl              wn         cm-1           # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[um or nm])
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30
```

```bash
# Data source #
Database                                ExoMol
Molecule                                CO2
Isotopologue                            12C-16O2
Dataset                                 UCL-4000
MolIsoID                                21


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     1
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               4
ChunkSize                               1000000


# Calculate oscillator strengths #
gf/f                                    gf
PlotOscillatorStrength(Y/N)             Y         
PlotOscillatorStrengthMethod            lin                       # Plot in linear (lin) or logarithm (log)
PlotOscillatorStrengthWnWl              wl         um             # Wavenumber (wn in unit cm-1) or wavelength (wl in unit[um or nm])
Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30
```
