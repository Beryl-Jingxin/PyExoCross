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

If `UncFilter(Y/N)` is yes, the value is the maximum uncertainty you require.

If `Threshold(Y/N)` is yes, the value is the minimum intensity you require.

If `QNsFilter(Y/N)` is yes, the spelling of the quantum number labels must be the same as `QNslabel`.

Doppler profile uses Doppler HWHM calculated by program, if you want to use Doppler profile, set `N` after `DopplerHWHM(Y/N)`. If you use Gaussian profile, please set `Y` after `DopplerHWHM(Y/N)`, your Doppler HWHM value will be used for calculateing Gaussian profile.

If you want a figure of corss sections, please set `Y` for `PlotCrossSection(Y/N)`.

* The definition file `.def` of ExoMol database (available at [exomol.com](https://www.exomol.com/)) provides the labels and formats of the quantum numbers for each species for reference.
* HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.

**Note**

You can define the quantum number column name by yourself, but please make sure it has letters without any spaces.
e.g. 'c1', 'c2', 'v1', 'v2', 'electronicState', 'electronic_state', '1v', '2v', 'M/E/C'.
Wrong format of the quantum number column nams: '1', '2', 'electronic state'.

*Example*

```bash
UncFilter(Y/N)                          N          0.001          # If Y, default value 0.001
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
Cutoff(Y/N)                             Y          100            # If Y, default value 25
QNsFilter(Y/N)                          N          par[+]   e/f[e]   v[0,1,2,3]  
DopplerHWHM(Y/N)                        Y          0.1            # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     N          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   N
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
Broadeners                              Air    Self    
Ratios                                  0.7    0.3    
```

```bash
Broadeners                              H2    He   
Ratios                                  0.9   0.1   
```

## Line profiles

Choose line profile from:

`Doppler`, `Gaussian`, `Lorentzian`, `SciPyVoigt`, `SciPyWofzVoigt`, `HumlicekVoigt`, `PseudoVoigt`, `PseudoKielkopfVoigt`, `PseudoOliveroVoigt`, `PseudoLiuLinVoigt`, `PseudoRoccoVoigt`, `BinnedDoppler`, `BinnedGaussian`, `BinnedLorentzian`, `BinnedVoigt`.

*Example*

```bash
# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f

# Calculate stick spectra or cross sections #
Temperature                             300
Range                                   0          30000
Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
UncFilter(Y/N)                          Y          0.001          # If Y, default value 0.001
Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30

# Calculate cross sections #
Pressure                                1
Npoints/BinSize                         Npoints    10001
Broadeners                              Default    
Ratios                                  1.0        
Profile                                 Gaussian        
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'
Cutoff(Y/N)                             Y          100            # If Y, default value 25
QNsFilter(Y/N)                          Y          par[+]   e/f[e]   v[0,1,2,3]  
DopplerHWHM(Y/N)                        Y          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     Y          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   Y
```

```bash
# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f

# Calculate stick spectra or cross sections #
Temperature                             1000
Range                                   1000       5000
Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
UncFilter(Y/N)                          No          0.001         # If Y, default value 0.001
Threshold(Y/N)                          NO          1e-30         # If Y, default value 1e-30

# Calculate cross sections #
Pressure                                0.1
Npoints/BinSize                         BinSize    0.1
Broadeners                              Air    Self    
Ratios                                  0.7    0.3     
Profile                                 SciPyVoigt        
Wavenumber(wn)/wavelength(wl)           wl                        # 'wn' or 'wl'
Cutoff(Y/N)                             N          100            # If Y, default value 25
QNsFilter(Y/N)                          N          
DopplerHWHM(Y/N)                        n          0.1            # Set Doppler HWHM as a constant 
LorentzianHWHM(Y/N)                     n          0.5            # Set Lorentzian HWHM as a constant
PlotCrossSection(Y/N)                   n
```

**Note**

If the line lists data is not in the ExoMol format, please convert your data format into the ExoMol format at first and then compute partition functions with *PyExoCross*.
