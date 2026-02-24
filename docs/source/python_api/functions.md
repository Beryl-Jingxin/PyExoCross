# Functions Reference

All functions are accessible from the top-level `pyexocross` namespace:

```python
import pyexocross as px
```

Every function accepts either:
- An `.inp` file path via `inp_filepath`, **or**
- Keyword arguments (`**kwargs`), **or**
- Both (keyword arguments override `.inp` file values).

---

## `px.run`

```python
px.run(inp_filepath, force_reload=False)
```

Run **all** enabled functions from an `.inp` configuration file.  This is the
programmatic equivalent of running `python3 run.py -p <inp_filepath>` from the
command line.

**Parameters**

| Parameter | Type | Required | Description |
|---|---|---|---|
| `inp_filepath` | `str` | Yes | Path to the `.inp` configuration file |
| `force_reload` | `bool` | No | Force re-parse of `.inp` even when cached in current Python process (default: `False`) |

:::{note}
If you edit the same `.inp` file repeatedly in a long-lived session (e.g. Jupyter),
use `force_reload=True` to avoid reusing stale cached parameters.
:::

**Example**

```python
px.run('/path/to/MgH_ExoMol.inp')

# Re-read updated .inp content in the same session
px.run('/path/to/MgH_ExoMol.inp', force_reload=True)
```

---

## `px.conversion`

```python
px.conversion(inp_filepath=None, **kwargs)
```

Convert between ExoMol/ExoAtom and HITRAN line-list formats.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inp_filepath` | `str` | `None` | Path to `.inp` file |
| `database` | `str` | `'ExoMol'` | `'ExoMol'`, `'ExoAtom'`, `'HITRAN'`, or `'HITEMP'` |
| `molecule` | `str` | `None` | Molecule name (e.g. `'MgH'`, `'NO'`) |
| `atom` | `str` | `None` | Atom name for ExoAtom (e.g. `'Ar'`, `'Li'`) |
| `isotopologue` | `str` | `None` | Isotopologue (e.g. `'24Mg-1H'`, `'14N-16O'`) |
| `dataset` | `str` | `None` | Dataset name (e.g. `'XAB'`, `'NIST'`) |
| `species_id` | `int` | `0` | Species identifier (e.g. `501`, `81`) |
| `read_path` | `str` | `'./'` | Path to input database |
| `save_path` | `str` | `'./output/'` | Path for output files |
| `logs_path` | `str` | `None` | Log file path |
| `conversion_format` | `int` | `1` | `1` = ExoMol/ExoAtom -> HITRAN/HITEMP, `2` = HITRAN/HITEMP -> ExoMol/ExoAtom |
| `conversion_min_freq` | `float` | `0` | Minimum frequency in cm⁻¹ |
| `conversion_max_freq` | `float` | `30000` | Maximum frequency in cm⁻¹ |
| `conversion_unc` | `float` or `None` | `None` | Uncertainty filter (cm⁻¹); `None` to disable |
| `conversion_threshold` | `float` or `None` | `None` | Intensity threshold (cm/molecule); `None` to disable |
| `global_qn_label_list` | `list[str]` | `[]` | Global quantum number labels for HITRAN output |
| `global_qn_format_list` | `list[str]` | `[]` | Global QN format specifiers (total 15 chars in HITRAN2004) |
| `local_qn_label_list` | `list[str]` | `[]` | Local quantum number labels for HITRAN output |
| `local_qn_format_list` | `list[str]` | `[]` | Local QN format specifiers (total 15 chars in HITRAN2004) |
| `qnslabel_list` | `list[str]` | `[]` | Quantum number labels from `.states` file |
| `qnsformat_list` | `list[str]` | `[]` | Quantum number format specifiers |
| `ncputrans` | `int` | `4` | CPU cores for processing transitions |
| `ncpufiles` | `int` | `1` | Files processed simultaneously |
| `chunk_size` | `int` | `100000` | Chunk size for transitions |

**Legacy Aliases**

```python
px.convert_exomol_to_hitran(...)  # Sets conversion_format=1
px.convert_hitran_to_exomol(...)  # Sets conversion_format=2
```

**Example**

```python
px.conversion(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    conversion_format=1,
    global_qn_label_list=['eS', 'v', 'Omega'],
    global_qn_format_list=['%9s', '%2d', '%4s'],
    local_qn_label_list=['J', 'e/f'],
    local_qn_format_list=['%5.1f', '%2s'],
)
```

---

## `px.partition_functions`

```python
px.partition_functions(inp_filepath=None, **kwargs)
```

Calculate partition functions $Q(T) = \sum_n g_n^{\text{tot}} \exp(-c_2 \tilde{E}_n / T)$.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inp_filepath` | `str` | `None` | Path to `.inp` file |
| `database` | `str` | `'ExoMol'` | Database name |
| `molecule` | `str` | `None` | Molecule name |
| `atom` | `str` | `None` | Atom name (ExoAtom) |
| `isotopologue` | `str` | `None` | Isotopologue |
| `dataset` | `str` | `None` | Dataset name |
| `species_id` | `int` | `0` | Species identifier |
| `read_path` | `str` | `'./'` | Input database path |
| `save_path` | `str` | `'./output/'` | Output path |
| `logs_path` | `str` | `None` | Log file path |
| `ntemp` | `int` | `1` | Temperature step in K |
| `tmax` | `int` | `5000` | Maximum temperature in K |
| `ncputrans` | `int` | `4` | CPU cores for transitions |
| `ncpufiles` | `int` | `1` | Files processed simultaneously |
| `chunk_size` | `int` | `100000` | Chunk size |

**Legacy Alias**: `px.partition_function`

**Example**

```python
px.partition_functions(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    ntemp=1,
    tmax=5000,
)
```

---

## `px.specific_heats`

```python
px.specific_heats(inp_filepath=None, **kwargs)
```

Calculate specific heat capacities $C_p(T)$ using partition function derivatives.

**Parameters**: Same as [`px.partition_functions`](#pxpartition_functions).

**Legacy Alias**: `px.specific_heat`

**Example**

```python
px.specific_heats(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    ntemp=1,
    tmax=5000,
)
```

---

## `px.cooling_functions`

```python
px.cooling_functions(inp_filepath=None, **kwargs)
```

Calculate cooling functions:

$$W(T) = \frac{1}{4\pi Q(T)} \sum_{f,i} A_{fi}\, h c \tilde{\nu}_{fi}\, g'\, \exp\!\left(-\frac{c_2 \tilde{E}'}{T}\right)$$

**Parameters**: Same as [`px.partition_functions`](#pxpartition_functions).

**Legacy Alias**: `px.cooling_function`

**Example**

```python
px.cooling_functions(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    ntemp=1,
    tmax=5000,
)
```

---

## `px.lifetimes`

```python
px.lifetimes(inp_filepath=None, **kwargs)
```

Calculate radiative lifetimes:
$$\tau_i = 1 / \sum_f A_{fi}$$

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inp_filepath` | `str` | `None` | Path to `.inp` file |
| `database` | `str` | `'ExoMol'` | Database name |
| `molecule` | `str` | `None` | Molecule name |
| `atom` | `str` | `None` | Atom name (ExoAtom) |
| `isotopologue` | `str` | `None` | Isotopologue |
| `dataset` | `str` | `None` | Dataset name |
| `species_id` | `int` | `0` | Species identifier |
| `read_path` | `str` | `'./'` | Input database path |
| `save_path` | `str` | `'./output/'` | Output path |
| `logs_path` | `str` | `None` | Log file path |
| `compress` | `bool` | `False` | `True` to save as `.bz2`; `False` for uncompressed |
| `ncputrans` | `int` | `4` | CPU cores for transitions |
| `ncpufiles` | `int` | `1` | Files processed simultaneously |
| `chunk_size` | `int` | `100000` | Chunk size |

**Legacy Alias**: `px.lifetime`

**Example**

```python
px.lifetimes(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    compress=False,
)
```

---

## `px.oscillator_strengths`

```python
px.oscillator_strengths(inp_filepath=None, **kwargs)
```

Calculate oscillator strengths: weighted $gf = g_{\text{tot}}' A_{fi} / (c \tilde{\nu}_{fi})^2$ or
unweighted $f$-values.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inp_filepath` | `str` | `None` | Path to `.inp` file |
| `database` | `str` | `'ExoMol'` | Database name |
| `molecule` | `str` | `None` | Molecule name |
| `atom` | `str` | `None` | Atom name (ExoAtom) |
| `isotopologue` | `str` | `None` | Isotopologue |
| `dataset` | `str` | `None` | Dataset name |
| `species_id` | `int` | `0` | Species identifier |
| `read_path` | `str` | `'./'` | Input database path |
| `save_path` | `str` | `'./output/'` | Output path |
| `logs_path` | `str` | `None` | Log file path |
| `gf_or_f` | `str` | `'f'` | `'gf'` for weighted, `'f'` for unweighted |
| `ncputrans` | `int` | `4` | CPU cores for transitions |
| `ncpufiles` | `int` | `1` | Files processed simultaneously |
| `chunk_size` | `int` | `100000` | Chunk size |
| `plot` | `bool` | `False` | Whether to generate a plot |
| `plot_method` | `str` | `'log'` | `'log'` or `'linear'` |
| `plot_wn_wl` | `str` | `'WN'` | `'WN'` (wavenumber) or `'WL'` (wavelength) |
| `plot_unit` | `str` | `'cm-1'` | `'cm-1'`, `'um'`, or `'nm'` |
| `limit_yaxis` | `float` | `1e-30` | Lower limit for y-axis |

**Legacy Alias**: `px.oscillator_strength`

**Example**

```python
px.oscillator_strengths(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    gf_or_f='f',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN',
    plot_unit='cm-1',
    limit_yaxis=1e-30,
)
```

---

## `px.stick_spectra`

```python
px.stick_spectra(inp_filepath=None, **kwargs)
```

Calculate LTE or Non-LTE stick spectra (absorption or emission).

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inp_filepath` | `str` | `None` | Path to `.inp` file |
| **Species** | | | |
| `database` | `str` | `'ExoMol'` | `'ExoMol'`, `'ExoAtom'`, `'HITRAN'`, or `'HITEMP'` |
| `molecule` | `str` | `None` | Molecule name |
| `atom` | `str` | `None` | Atom name (ExoAtom) |
| `isotopologue` | `str` | `None` | Isotopologue |
| `dataset` | `str` | `None` | Dataset name |
| `species_id` | `int` | `0` | Species identifier |
| **Paths** | | | |
| `read_path` | `str` | `'./'` | Input database path |
| `save_path` | `str` | `'./output/'` | Output path |
| `logs_path` | `str` | `None` | Log file path |
| **Physical conditions** | | | |
| `temperatures` | `list[float]` | `[1000]` | Temperature(s) in K |
| `wn_wl` | `str` | `'WN'` | `'WN'` (wavenumber) or `'WL'` (wavelength) |
| `wn_wl_unit` | `str` | `'cm-1'` | `'cm-1'`, `'um'`, or `'nm'` |
| `min_range` | `float` | `0` | Minimum range value (in `wn_wl_unit`) |
| `max_range` | `float` | `30000` | Maximum range value (in `wn_wl_unit`) |
| `abs_emi` | `str` | `'Ab'` | `'Ab'` for absorption, `'Em'` for emission |
| **Non-LTE** | | | |
| `nlte_method` | `str` | `'L'` | `'L'`=LTE, `'T'`=Treanor, `'D'`=Density, `'P'`=Population |
| `tvib_list` | `list[float]` | `[]` | Vibrational temperatures for method `'T'` |
| `trot_list` | `list[float]` | `[]` | Rotational temperatures for method `'T'` or `'D'` |
| `vib_label` | `list[str]` | `[]` | Vibrational QN labels for method `'T'` |
| `rot_label` | `list[str]` | `[]` | Rotational QN labels for method `'T'` |
| `nlte_path` | `str` | `None` | Path to NLTE data file (methods `'D'` and `'P'`) |
| **Filters** | | | |
| `threshold` | `float` or `None` | `None` | Intensity threshold (cm/molecule); `None` to disable |
| `unc_filter` | `float` or `None` | `None` | Uncertainty filter (cm⁻¹); `None` to disable |
| `qns_filter` | `dict` or `None` | `None` | Quantum number filter (see below) |
| **Quantum numbers** | | | |
| `qnslabel_list` | `list[str]` | `[]` | QN labels from `.states` file |
| `qnsformat_list` | `list[str]` | `[]` | QN format specifiers |
| **Computing** | | | |
| `ncputrans` | `int` | `4` | CPU cores for transitions |
| `ncpufiles` | `int` | `1` | Files processed simultaneously |
| `chunk_size` | `int` | `100000` | Chunk size |
| **Plotting** | | | |
| `plot` | `bool` | `False` | Whether to generate a plot |
| `plot_method` | `str` | `'log'` | `'log'` or `'linear'` |
| `plot_wn_wl` | `str` | `'WN'` | X-axis: `'WN'` or `'WL'` |
| `plot_unit` | `str` | `'cm-1'` | `'cm-1'`, `'um'`, or `'nm'` |
| `limit_yaxis` | `float` | `1e-30` | Lower limit for y-axis (cm/molecule) |

### Quantum Number Filter (`qns_filter`)

The `qns_filter` parameter accepts a dictionary where keys are quantum number
labels and values are lists of accepted value patterns.  Use an empty list
`[]` to accept all values for that label.

```python
qns_filter={
    'par': [],                    # accept all parities
    'v': ['0,', '1,', '2,'],     # only v = 0, 1, 2
    'eS': [],                    # accept all electronic states
}
```

**Example**

```python
px.stick_spectra(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    temperatures=[1000, 2000],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    abs_emi='Ab',
)
```

---

## `px.cross_sections`

```python
px.cross_sections(inp_filepath=None, **kwargs)
```

Calculate LTE or Non-LTE cross sections with a choice of line profiles.

**Parameters**

All parameters from [`px.stick_spectra`](#pxstick_spectra), plus:

| Parameter | Type | Default | Description |
|---|---|---|---|
| **Grid** | | | |
| `pressures` | `list[float]` | `[1.0]` | Pressure(s) in bar |
| `bin_size` | `float` | `0.1` | Bin size (in `wn_wl_unit`); mutually exclusive with `n_point` |
| `n_point` | `int` | `None` | Number of grid points; mutually exclusive with `bin_size` |
| **Line profile** | | | |
| `profile` | `str` | `'Gaussian'` | Line profile (see table below) |
| `cutoff` | `float` or `None` | `None` | Wing cutoff in cm⁻¹; `None` to disable |
| `predissociation` | `bool` | `False` | Use predissociation lifetimes for Voigt width |
| **Broadening** | | | |
| `broadeners` | `list[str]` | `['Default']` | Broadening species |
| `ratios` | `list[float]` | `[1.0]` | Broadener mixing ratios (must sum to 1.0) |
| `alpha_hwhm` | `float` or `None` | `3.0` | Constant Doppler HWHM; `None` to calculate |
| `gamma_hwhm` | `float` or `None` | `None` | Constant Lorentzian HWHM; `None` to calculate |
| **Plotting** | | | |
| `plot` | `bool` | `False` | Whether to generate a plot |
| `plot_method` | `str` | `'log'` | `'log'` or `'linear'` |
| `plot_wn_wl` | `str` | `'WN'` | X-axis: `'WN'` or `'WL'` |
| `plot_unit` | `str` | `'cm-1'` | `'cm-1'`, `'um'`, or `'nm'` |
| `limit_yaxis` | `float` | `1e-30` | Lower limit for y-axis (cm$^2$/molecule) |

### Available Line Profiles

| Profile Name | Description |
|---|---|
| `'Doppler'` | Pure Doppler profile |
| `'Gaussian'` | Gaussian profile |
| `'Lorentzian'` | Lorentzian profile |
| `'SciPyVoigt'` | Voigt via SciPy `wofz` (recommended) |
| `'SciPyWofzVoigt'` | Voigt via SciPy `wofz` (alternative) |
| `'HumlicekVoigt'` | Humlicek algorithm for Voigt |
| `'ThompsonPseudoVoigt'` | Thompson pseudo-Voigt approximation |
| `'KielkopfPseudoVoigt'` | Kielkopf pseudo-Voigt approximation |
| `'OliveroPseudoVoigt'` | Olivero pseudo-Voigt approximation |
| `'LiuLinPseudoVoigt'` | Liu-Lin pseudo-Voigt approximation |
| `'RoccoPseudoVoigt'` | Rocco pseudo-Voigt approximation |
| `'BinnedDoppler'` | Binned Doppler profile |
| `'BinnedGaussian'` | Binned Gaussian profile |
| `'BinnedLorentzian'` | Binned Lorentzian profile |
| `'BinnedVoigt'` | Binned Voigt profile |

### Available Broadeners

| Broadener | Description |
|---|---|
| `'Default'` | Default broadening parameters |
| `'Air'` | Air broadening (N$_2$ + O$_2$) |
| `'Self'` | Self broadening |
| `'H2'` | Hydrogen broadening |
| `'He'` | Helium broadening |
| `'CO2'` | CO$_2$ broadening |
| `'H2O'` | Water broadening |

**Legacy Alias**: `px.cross_section`

**Example**

```python
px.cross_sections(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    temperatures=[1000, 2000],
    pressures=[1.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    bin_size=0.1,
    profile='SciPyVoigt',
    broadeners=['Default'],
    ratios=[1.0],
    cutoff=25.0,
    plot=True,
    plot_method='log',
)
```
