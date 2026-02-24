# Parameters Reference

This page provides a **complete** reference for all keyword arguments accepted
by PyExoCross API functions.  Parameters are organized by category.

---

## Species Identification

| Parameter | Type | Default | Description | Databases |
|---|---|---|---|---|
| `database` | `str` | `'ExoMol'` | Database type | All |
| `molecule` | `str` | `None` | Molecule name (e.g. `'H2O'`, `'MgH'`, `'NO'`) | ExoMol, HITRAN, HITEMP |
| `atom` | `str` | `None` | Atom name (e.g. `'Ar'`, `'Li'`, `'Al'`) | ExoAtom |
| `isotopologue` | `str` | `None` | Isotopologue (e.g. `'1H2-16O'`, `'24Mg-1H'`) | ExoMol, HITRAN, HITEMP |
| `dataset` | `str` | `None` | Dataset name (e.g. `'POKAZATEL'`, `'XAB'`, `'NIST'`) | All |
| `species_id` | `int` | `0` | Numeric species identifier | All |

### `database` Accepted Values

| Value | Description |
|---|---|
| `'ExoMol'` | ExoMol molecular database (default) |
| `'ExoAtom'` | ExoAtom atomic database |
| `'HITRAN'` | HITRAN molecular spectroscopic database |
| `'HITEMP'` | HITEMP high-temperature extension of HITRAN |

### `species_id` Encoding

The `species_id` is a compact numeric identifier.  The first digit(s) encode
the species main ID and the last digit encodes the species sub ID:

```
species_id = species_main_id * 10 + species_sub_id
```

| Database | Example | `species_id` | Main ID | Sub ID | Meaning |
|---|---|---|---|---|---|
| ExoMol | MgH (24Mg-1H) | `501` | `50` | `1` | HITRAN molecule ID #50, isotopologue ID #1 |
| ExoAtom | Ar | `601` | `60` | `1` | HITRAN atom ID #60, isotope ID #1 |
| HITRAN | H2O (1H2-16O) | `11` | `1` | `1` | HITRAN molecule ID #1, isotopologue ID #1 |
| HITRAN | CO2 (12C-16O2) | `21` | `2` | `1` | HITRAN molecule ID #2, isotopologue ID #1 |
| HITRAN | NO (14N-16O) | `81` | `8` | `1` | HITRAN molecule ID #8, isotopologue ID #1 |
| HITRAN | C2H2 | `261` | `26` | `1` | HITRAN molecule ID #26, isotopologue ID #1 |

---

## File Paths

| Parameter | Type | Default | Description |
|---|---|---|---|
| `read_path` | `str` | `'./'` | Path to input data |
| `save_path` | `str` | `'./output/'` | Path for output files |
| `logs_path` | `str` | `None` | Log file path; `None` to skip logging |

### `read_path` Conventions

| Database | `read_path` Format | Example |
|---|---|---|
| ExoMol  | Root directory of ExoMol database  | `'/data/ExoMol/'`  |
| ExoAtom | Root directory of ExoAtom database | `'/data/ExoAtom/'` |
| HITRAN  | Root directory of HITRAN database  | `'/data/HITRAN/'`  |
| HITEMP  | Root directory of HITEMP database  | `'/data/HITEMP/'`  |

For ExoMol, the code automatically resolves the full path as:

```
{read_path}/{molecule}/{isotopologue}/{dataset}/
```

For ExoAtom:

```
{read_path}/{atom}/{dataset}/
```

For HITRAN/HITEMP:

```
  {read_path}/{molecule}/{isotopologue}/
```

---

## Computing Resources

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ncputrans`  | `int` | `4` | Number of CPU cores for processing each transitions file |
| `ncpufiles`  | `int` | `1` | Number of transitions files processed simultaneously |
| `chunk_size` | `int` | `100000` | Chunk size when reading/calculating transitions |

:::{tip}
For large line lists (e.g. H2O with billions of lines), increase `chunk_size`
to reduce I/O overhead.  Set `ncpufiles > 1` if the species has many
separate `.trans` files (common for ExoMol).
:::

---

## Conversion Parameters

Used by `px.conversion()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `conversion_format` | `int` | `1` | `1` = ExoMol/ExoAtom -> HITRAN, `2` = HITRAN -> ExoMol |
| `conversion_min_freq` | `float` | `0` | Minimum frequency in cm⁻¹ |
| `conversion_max_freq` | `float` | `30000` | Maximum frequency in cm⁻¹ |
| `conversion_unc` | `float` or `None` | `None` | Max uncertainty (cm⁻¹); `None` disables |
| `conversion_threshold` | `float` or `None` | `None` | Min intensity (cm/molecule); `None` disables |
| `global_qn_label_list` | `list[str]` | `[]` | Global QN labels for HITRAN output |
| `global_qn_format_list` | `list[str]` | `[]` | Global QN Fortran-style formats |
| `local_qn_label_list` | `list[str]` | `[]` | Local QN labels for HITRAN output |
| `local_qn_format_list` | `list[str]` | `[]` | Local QN Fortran-style formats |

### Quantum Number Formats

The global and local quantum number format strings follow HITRAN2004 convention.
Each format field is a C/Fortran-style specifier:

```python
global_qn_label_list  = ['eS',  'v',  'Omega']
global_qn_format_list = ['%9s', '%2d', '%4s']    # Total: 15 chars

local_qn_label_list   = ['J',    'e/f']
local_qn_format_list  = ['%5.1f', '%2s']         # Total: 15 chars (padded)
```

:::{important}
For HITRAN2004 format, both global and local quantum number fields are
**exactly 15 characters** each.  Ensure the sum of format widths matches.
:::

---

## Partition Function / Specific Heat / Cooling Function Parameters

Used by `px.partition_functions()`, `px.specific_heats()`, `px.cooling_functions()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `ntemp` | `int` | `1` | Temperature step interval in K |
| `tmax` | `int` | `5000` | Maximum temperature in K |

The calculation produces values at every `ntemp` K from 1 K to `tmax` K.

---

## Lifetime Parameters

Used by `px.lifetimes()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `compress` | `bool` | `False` | `True` saves output as `.states.bz2`; `False` as `.states` |

---

## Oscillator Strength Parameters

Used by `px.oscillator_strengths()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `gf_or_f` | `str` | `'f'` | `'gf'` for weighted oscillator strength, `'f'` for $f$-value |

---

## Physical Conditions

Used by `px.stick_spectra()` and `px.cross_sections()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `temperatures` | `list[float]` | `[1000]` | Temperature(s) in K |
| `pressures` | `list[float]` | `[1.0]` | Pressure(s) in bar (cross sections only) |
| `wn_wl` | `str` | `'WN'` | `'WN'` = wavenumber, `'WL'` = wavelength |
| `wn_wl_unit` | `str` | `'cm-1'` | Unit: `'cm-1'`, `'um'` (micron), or `'nm'` |
| `min_range` | `float` | `0` | Minimum of range (in `wn_wl_unit`) |
| `max_range` | `float` | `30000` | Maximum of range (in `wn_wl_unit`) |
| `abs_emi` | `str` | `'Ab'` | `'Ab'` = absorption, `'Em'` = emission |

---

## Non-LTE Parameters

Used by `px.stick_spectra()` and `px.cross_sections()` when `nlte_method != 'L'`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `nlte_method` | `str` | `'L'` | Non-LTE method (see table below) |
| `tvib_list` | `list[float]` | `[]` | Vibrational temperature(s) in K |
| `trot_list` | `list[float]` | `[]` | Rotational temperature(s) in K |
| `vib_label` | `list[str]` | `[]` | Vibrational quantum number labels |
| `rot_label` | `list[str]` | `[]` | Rotational quantum number labels |
| `nlte_path` | `str` | `None` | Path to NLTE data file |

### Non-LTE Methods

| Value | Name | Description | Required Parameters |
|---|---|---|---|
| `'L'` | LTE | Local thermodynamic equilibrium (default) | `temperatures` |
| `'T'` | Treanor | Two-temperature model | `temperatures`, `tvib_list`, `trot_list`, `vib_label`, `rot_label` |
| `'D'` | Density | Number density NLTE | `temperatures`, `trot_list`, `nlte_path` |
| `'P'` | Population | Population NLTE | `temperatures`, `nlte_path` |

---

## Filter Parameters

Used by `px.stick_spectra()` and `px.cross_sections()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `threshold` | `float` or `None` | `None` | Minimum intensity threshold (cm/molecule) |
| `unc_filter` | `float` or `None` | `None` | Maximum uncertainty filter (cm⁻¹) |
| `qns_filter` | `dict` or `None` | `None` | Quantum number filter dictionary |

### `qns_filter` Format

A dictionary mapping quantum number label to a list of accepted value patterns.
An empty list `[]` means "accept all values".

```python
# Accept only v = 0, 1, 2 and all parities
qns_filter = {
    'v': ['0,', '1,', '2,'],
    'par': [],
}

# Accept all (equivalent to None)
qns_filter = {
    'configuration': [],
    'Multiple': [],
    'parity': [],
}
```

---

## Quantum Number Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `qnslabel_list` | `list[str]` | `[]` | Quantum number column labels |
| `qnsformat_list` | `list[str]` | `[]` | Quantum number format specifiers |

The labels and formats must correspond to the columns in the `.states` file
**after** the standard columns (ID, energy, degeneracy, lifetime).

**ExoMol Example** (MgH XAB):

```python
qnslabel_list  = ['par', 'e/f', 'eS', 'v', 'Lambda', 'Sigma', 'Omega']
qnsformat_list = ['%1s', '%1s', '%13s', '%3d', '%2d', '%7.1f', '%7.1f']
```

**ExoAtom Example** (Ar NIST):

```python
qnslabel_list  = ['configuration', 'Multiple', 'parity']
qnsformat_list = ['%20s', '%10s', '%2s']
```

**HITRAN Example** (NO):

```python
qnslabel_list  = ['J', 'X', 'Omega', 'v1', 'Sym', 'F']
qnsformat_list = ['%5.1f', '%2s', '%3s', '%2d', '%1s', '%5s']
```

---

## Line Profile Parameters

Used by `px.cross_sections()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `profile` | `str` | `'Gaussian'` | Line profile type |
| `cutoff` | `float` or `None` | `None` | Wing cutoff distance (cm⁻¹); `None` for no cutoff |
| `predissociation` | `bool` | `False` | Include predissociation lifetimes for Voigt width |

---

## Broadening Parameters

Used by `px.cross_sections()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `broadeners` | `list[str]` | `['Default']` | Broadening species (e.g. `['H2', 'He']`) |
| `ratios` | `list[float]` | `[1.0]` | Mixing ratios (must sum to 1.0) |
| `alpha_hwhm` | `float` or `None` | `3.0` | Constant Doppler HWHM (cm⁻¹); `None` to auto-calculate |
| `gamma_hwhm` | `float` or `None` | `None` | Constant Lorentzian HWHM (cm⁻¹); `None` to auto-calculate |

:::{tip}
Set `alpha_hwhm=None` and `gamma_hwhm=None` to let PyExoCross calculate
HWHMs from broadening parameters in the database.  Use constant values for
quick testing or when database broadening files are unavailable.
:::

---

## Grid Parameters

Used by `px.cross_sections()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `bin_size` | `float` | `0.1` | Bin size in `wn_wl_unit` (mutually exclusive with `n_point`) |
| `n_point` | `int` | `None` | Number of grid points (mutually exclusive with `bin_size`) |

:::{note}
If both `bin_size` and `n_point` are specified, `n_point` takes precedence
and `bin_size` is calculated from the range.
:::

---

## Plotting Parameters

Used by `px.oscillator_strengths()`, `px.stick_spectra()`, `px.cross_sections()`.

| Parameter | Type | Default | Description |
|---|---|---|---|
| `plot` | `bool` | `False` | Whether to generate a plot |
| `plot_method` | `str` | `'log'` | `'log'` or `'linear'` scale |
| `plot_wn_wl` | `str` | `'WN'` | X-axis: `'WN'` (wavenumber) or `'WL'` (wavelength) |
| `plot_unit` | `str` | `'cm-1'` | X-axis unit: `'cm-1'`, `'um'`, or `'nm'` |
| `limit_yaxis` | `float` | `1e-30` | Lower limit for y-axis |

:::{note}
These are **convenience** kwargs that get remapped internally to
function-specific names.  For example, `plot=True` in
`px.cross_sections()` becomes `plot_cross_section=True`.
:::
