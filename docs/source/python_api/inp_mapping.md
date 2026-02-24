# `.inp` File to Python Kwarg Mapping

This page provides a complete mapping between the `.inp` configuration file
keywords and the corresponding Python keyword arguments.  Use this as a
reference when migrating from `.inp`-based workflows to the Python API.

---

## Species Identification

| `.inp` Keyword | Python Kwarg | Example (`.inp`) | Example (Python) |
|---|---|---|---|
| `Database` | `database` | `ExoMol` | `'ExoMol'` |
| `Molecule` | `molecule` | `MgH` | `'MgH'` |
| `Atom` | `atom` | `Ar` | `'Ar'` |
| `Isotopologue` | `isotopologue` | `24Mg-1H` | `'24Mg-1H'` |
| `Dataset` | `dataset` | `XAB` | `'XAB'` |
| `SpeciesID` | `species_id` | `501` | `501` |

---

## File Paths

| `.inp` Keyword | Python Kwarg | Example (`.inp`) | Example (Python) |
|---|---|---|---|
| `ReadPath` | `read_path` | `/data/ExoMol/` | `'/data/ExoMol/'` |
| `SavePath` | `save_path` | `/data/output/` | `'/data/output/'` |
| `LogFilePath` | `logs_path` | `/data/output/log/run.log` | `'/data/output/log/run.log'` |

---

## Function Toggles

| `.inp` Keyword | Python API Function | `.inp` Value | Python Usage |
|---|---|---|---|
| `Conversion` | `px.conversion()` | `1` | Call function directly |
| `PartitionFunctions` | `px.partition_functions()` | `1` | Call function directly |
| `SpecificHeats` | `px.specific_heats()` | `1` | Call function directly |
| `CoolingFunctions` | `px.cooling_functions()` | `1` | Call function directly |
| `Lifetimes` | `px.lifetimes()` | `1` | Call function directly |
| `OscillatorStrengths` | `px.oscillator_strengths()` | `1` | Call function directly |
| `StickSpectra` | `px.stick_spectra()` | `1` | Call function directly |
| `CrossSections` | `px.cross_sections()` | `1` | Call function directly |

:::{note}
In the Python API, you don't need to set these toggle flags.  Simply calling
a function (e.g. `px.cross_sections(...)`) automatically enables that
function.  The toggle is only relevant in `.inp` files or when using
`px.run()`.
:::

---

## Computing Resources

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `NCPUtrans` | `ncputrans` | `2` | `2` (int) |
| `NCPUfiles` | `ncpufiles` | `4` | `4` (int) |
| `ChunkSize` | `chunk_size` | `1000000` | `1000000` (int) |

---

## Quantum Numbers

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `QNslabel` | `qnslabel_list` | `par e/f eS v Lambda Sigma Omega` | `['par', 'e/f', 'eS', 'v', 'Lambda', 'Sigma', 'Omega']` |
| `QNsformat` | `qnsformat_list` | `%1s %1s %13s %3d %2d %7.1f %7.1f` | `['%1s', '%1s', '%13s', '%3d', '%2d', '%7.1f', '%7.1f']` |

---

## Conversion

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `ConversionFormat` | `conversion_format` | `1` | `1` (int) |
| `ConversionFrequncyRange` | `conversion_min_freq`, `conversion_max_freq` | `0 30000` | Two separate `float` args |
| `GlobalQNLabel` | `global_qn_label_list` | `eS v Omega` | `['eS', 'v', 'Omega']` |
| `GlobalQNFormat` | `global_qn_format_list` | `%9s %2d %4s` | `['%9s', '%2d', '%4s']` |
| `LocalQNLabel` | `local_qn_label_list` | `J e/f` | `['J', 'e/f']` |
| `LocalQNFormat` | `local_qn_format_list` | `%5.1f %2s` | `['%5.1f', '%2s']` |
| `ConvUncFilter(Y/N)` | `conversion_unc` | `Y 0.01` | `0.01` (or `None` for `N`) |
| `ConvThreshold(Y/N)` | `conversion_threshold` | `Y 1e-30` | `1e-30` (or `None` for `N`) |

---

## Partition Functions / Specific Heats / Cooling Functions

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `Ntemp` | `ntemp` | `1` | `1` (int) |
| `Tmax` | `tmax` | `5000` | `5000` (int) |

---

## Lifetimes

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `Compress(Y/N)` | `compress` | `Y` or `N` | `True` or `False` |

---

## Oscillator Strengths

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `gf/f` | `gf_or_f` | `gf` or `f` | `'gf'` or `'f'` |
| `PlotOscillatorStrength(Y/N)` | `plot` | `Y` or `N` | `True` or `False` |
| `PlotOscillatorStrengthMethod` | `plot_method` | `log` | `'log'` or `'linear'` |
| `PlotOscillatorStrengthWnWl` | `plot_wn_wl`, `plot_unit` | `wn cm-1` | `'WN'`, `'cm-1'` |
| `Y-axisLimitOscillatorStrength` | `limit_yaxis` | `1e-30` | `float` |

---

## Stick Spectra

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `LTE/Non-LTE` | `nlte_method` | `LTE` or `Non-LTE` | `'L'` (LTE) or `'T'`/`'D'`/`'P'` |
| `Temperatures` | `temperatures` | `1000` | `[1000]` (list) |
| `WnWlUnit` | `wn_wl`, `wn_wl_unit` | `wn cm-1` | `'WN'`, `'cm-1'` |
| `Range` | `min_range`, `max_range` | `0 30000` | Two separate `float` args |
| `Absorption/Emission` | `abs_emi` | `Absorption` | `'Ab'` or `'Em'` |
| `UncFilter(Y/N)` | `unc_filter` | `Y 0.01` | `0.01` (or `None` for `N`) |
| `Threshold(Y/N)` | `threshold` | `Y 1e-30` | `1e-30` (or `None` for `N`) |
| `QNsFilter(Y/N)` | `qns_filter` | See `.inp` docs | `dict` (or `None`) |
| `PlotStickSpectra(Y/N)` | `plot` | `Y` | `True` or `False` |
| `PlotStickSpectraMethod` | `plot_method` | `log` | `'log'` or `'linear'` |
| `PlotStickSpectraWnWl` | `plot_wn_wl`, `plot_unit` | `wn cm-1` | `'WN'`, `'cm-1'` |
| `Y-axisLimitStickSpectra` | `limit_yaxis` | `1e-30` | `float` |

---

## Cross Sections

All stick spectra parameters above, plus:

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `Pressures` | `pressures` | `1.0` | `[1.0]` (list) |
| `Npoints/BinSize` | `n_point` or `bin_size` | `BinSize 0.1` | `bin_size=0.1` or `n_point=10000` |
| `Profile` | `profile` | `SciPyVoigt` | `'SciPyVoigt'` |
| `Cutoff(Y/N)` | `cutoff` | `Y 25` | `25.0` (or `None` for `N`) |
| `PredissocXsec(Y/N)` | `predissociation` | `Y` or `N` | `True` or `False` |
| `Broadeners` | `broadeners` | `H2 He` | `['H2', 'He']` |
| `Ratios` | `ratios` | `0.85 0.15` | `[0.85, 0.15]` |
| `DopplerHWHM(Y/N)` | `alpha_hwhm` | `Y 3.0` | `3.0` (or `None` for `N`) |
| `LorentzianHWHM(Y/N)` | `gamma_hwhm` | `Y 0.5` | `0.5` (or `None` for `N`) |
| `PlotCrossSection(Y/N)` | `plot` | `Y` | `True` or `False` |
| `PlotCrossSectionMethod` | `plot_method` | `log` | `'log'` or `'linear'` |
| `PlotCrossSectionWnWl` | `plot_wn_wl`, `plot_unit` | `wn cm-1` | `'WN'`, `'cm-1'` |
| `Y-axisLimitXsec` | `limit_yaxis` | `1e-30` | `float` |

---

## Non-LTE

| `.inp` Keyword | Python Kwarg | `.inp` Format | Python Type |
|---|---|---|---|
| `NLTEMethod` | `nlte_method` | `T` | `'T'` |
| `Tvib` | `tvib_list` | `1000 2000` | `[1000, 2000]` |
| `Trot` | `trot_list` | `300 400` | `[300, 400]` |
| `QNsVibLabel` | `vib_label` | `v, eS` | `['v', 'eS']` |
| `QNsRotLabel` | `rot_label` | `J, e/f` | `['J', 'e/f']` |
| *(NLTEMethod D/P path)* | `nlte_path` | `/path/to/file.csv` | `'/path/to/file.csv'` |

---

## Key Differences Between `.inp` and Python API

| Feature | `.inp` File | Python API |
|---|---|---|
| **Enable function** | Set flag to `1` | Call function directly |
| **Boolean values** | `Y` / `N` | `True` / `False` |
| **Disable filter** | `N` | `None` |
| **Lists** | Space-separated | Python `list` |
| **Multiple temperatures** | One per line | Single `list` argument |
| **Multiple pressures** | One per line | Single `list` argument |
| **Range** | Single line `min max` | Two args: `min_range`, `max_range` |
| **Plot toggle** | Function-specific name | Generic `plot=True` |
| **Unit spec** | Combined `wn cm-1` | Separate `wn_wl`, `wn_wl_unit` |

### Example Side-by-Side

**`.inp` file:**

```
Database                         ExoMol
Molecule                         MgH
Isotopologue                     24Mg-1H
Dataset                          XAB
SpeciesID                        501
CrossSections                    1
Temperatures                     1000:4000:1000
Pressures                        0.1,1.0
Range                            0  30000
Profile                          SciPyVoigt
```

**Python API equivalent:**

```python
px.cross_sections(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    temperatures=[1000, 2000, 3000, 4000],
    pressures=[0.1, 1.0],
    min_range=0,
    max_range=30000,
    profile='SciPyVoigt',
)
```
