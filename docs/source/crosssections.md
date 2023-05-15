Cross sections
==============

`Range`: Give two values as the minimum and maximum of the wavenumber range. No `,` or `;` between these two numbers, just leave blank here.

`Npoints/BinSize`: `Npoints` is the number of the points in grid. `BinSiza` is the interval size of the grid.

`Absorption/Emission`: Choose `Absorption` or `Emission`.

`Wavenumber(wn)/wavelength(wl)`: Choose `wn` or `wl`.

*Example*

```
Temperature                             300
Pressure                                1
Range                                   0          30000
Npoints/BinSize                         Npoints    30001

Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
```

## Filters

`Y/N`: `Y`, `YES`, `Yes`, `yes` and `N`, `NO`, `No`, `no` all can work. If you don't use it, write `N` here. You don't need to change the content behind it.

* The definition file `.def` of ExoMol database (available at [exomol.com](https://www.exomol.com/)) provides the labels and formats of the quantum numbers for each species for reference.
* HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.

**Note**

You can define the quantum number column name by yourself, but please make sure it has letters without any spaces.
e.g. 'c1', 'c2', 'v1', 'v2', 'electronicState', 'electronic_state', '1v', '2v', 'M/E/C'.
Wrong format of the quantum number column nams: '1', '2', 'electronic state'.

*Example*

```bash
Cutoff(Y/N)                             Y          25             # Default value 25
Threshold(Y/N)                          Y          1e-30          # Default value 1e-30
UncFilter(Y/N)                          Y          0.001
QNsFilter(Y/N)                          YES        par[+]   e/f[e]   v[0,1,2,3]
DopplerHWHM(Y/N)                        Yes        3              # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     yes        0.7            # Set Lorentzian HWHM as a constant
```

## Broadeners

`Default` means:

1. Default value of temperature exponent for all lines is 0.5.
2. Default value of Lorentzian half-width for all lines is 0.07.

If you want to use the default values of the broadeners, please set the `Ratios` of `Broadeners` under `Default` as `1.0`.

The broadening types and ratio values are corresponding, please do not change the place or delete elements.

*Example*

```bash
Broadeners                              Default    Air    Self    H2    He    CO2
Ratios                                  1.0        0.0    0.0     0.0   0.0   0.0
```

```bash
Broadeners                              Default    Air    Self    H2    He    CO2
Ratios                                  0.0        0.7    0.3     0.0   0.0   0.0
```

```bash
Broadeners                              Default    Air    Self    H2    He    CO2
Ratios                                  0.0        0.0    0.0     0.9   0.1   0.0
```

## Line profiles

Choose line profile from:

`Doppler`, `Gaussian`, `Lorentzian`, `SciPyVoigt`, `SciPyWofzVoigt`, `PseudoVoigt`, `PseudoKielkopfVoigt`, `PseudoOliveroVoigt`, `PseudoLiuLinVoigt`, `PseudoRoccoVoigt`, `BinnedDoppler`, `BinnedGaussian`, `BinnedLorentzian`, `BinnedVoigt`.

*Example*

```bash
# Calculate stick spectra or cross-sections #
Temperature                             300
Pressure                                1
Range                                   0          30000
Npoints/BinSize                         Npoints    30001

Cutoff(Y/N)                             Y          25             # Default value 25
Threshold(Y/N)                          Y          1e-30          # Default value 1e-30
UncFilter(Y/N)                          Y          0.001
QNsFilter(Y/N)                          YES        par[+]   e/f[e]   v[0,1,2,3]
DopplerHWHM(Y/N)                        Yes        3              # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     yes        0.7            # Set Lorentzian HWHM as a constant

Broadeners                              Default    Air    Self    H2    He    CO2
Ratios                                  1.0        0.0    0.0     0.0   0.0   0.0

Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
Profile                                 Gaussian
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
```

```bash
# Calculate stick spectra or cross-sections #
Temperature                             1000
Pressure                                1.2
Range                                   1000       5000
Npoints/BinSize                         BinSize    0.1

Cutoff(Y/N)                             N          25             # Default value 25
Threshold(Y/N)                          N          1e-30          # Default value 1e-30
UncFilter(Y/N)                          N          0.001
QNsFilter(Y/N)                          NO         par[+]   e/f[e]   v[0,1,2,3]
DopplerHWHM(Y/N)                        No         3              # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     no         0.7            # Set Lorentzian HWHM as a constant

Broadeners                              Default    Air    Self    H2    He    CO2
Ratios                                  0.0        0.7    0.3     0.0   0.0   0.0

Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
Profile                                 SciPyVoigt
Wavenumber(wn)/wavelength(wl)           wl                        # 'wn' or 'wl'
```

**Note**

If the line lists data is not in the ExoMol format, please convert your data format into the ExoMol format at first and then compute partition functions with *PyExoCross*.
