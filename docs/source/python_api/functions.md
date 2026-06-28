# Functions Reference

All functions are accessible from the top-level `pyexocross` namespace:

```python
import pyexocross as px
```

Every function accepts either:
- An `.inp` file path via `inp_filepath`, **or**
- Keyword arguments (`**kwargs`), **or**
- Both (keyword arguments override `.inp` file values).

## Common CPU/GPU Compute Kwargs

PyExoCross uses CPU by default.  You can explicitly select compute mode with:

| Kwarg | Type | Default | Description |
|---|---|---|---|
| `device` or `run_mode` | `str` | `'CPU'` | `'CPU'` or `'GPU'` |
| `gpu_backend` | `str` | `'AUTO'` | GPU backend policy: `'AUTO'`, `'CUDA'`, `'PyTorch-CUDA'`, `'CuPy-CUDA'`, or `'MPS'` |
| `gpu_batch_lines` | `int` | `8192` | GPU line-batch size for memory control |
| `gpu_batch_grid` | `int` | `256` | GPU grid-batch size for memory control |

GPU acceleration is available for:
- `px.cooling_functions` -> Recommend to use CPU mode
- `px.stick_spectra` -> Recommend to use CPU mode
- `px.cross_sections` 
- `px.stick_spectra_cross_section`

CPU formulas are used for:
- `px.conversion`
- `px.partition_functions`
- `px.specific_heats`
- `px.lifetimes`
- `px.oscillator_strengths`

If `device='GPU'` but no compatible backend is available, PyExoCross falls
back to CPU formulas.

```python
# Default CPU mode
px.cross_sections(..., device='CPU')

# GPU mode (auto backend)
px.cross_sections(
    ...,
    device='GPU',
    gpu_backend='AUTO',
    gpu_batch_lines=8192,
    gpu_batch_grid=256,
)

# CUDA policy: PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback
px.cross_sections(..., device='GPU', gpu_backend='CUDA')

# Force PyTorch CUDA only
px.cross_sections(..., device='GPU', gpu_backend='PyTorch-CUDA')

# Force CuPy CUDA only
px.cross_sections(..., device='GPU', gpu_backend='CuPy-CUDA')

# Force PyTorch MPS only
px.cross_sections(..., device='GPU', gpu_backend='MPS')
```

---

## `px.load`

```python
data = px.load(..., cache='auto', cache_dir=None, max_memory=512,
               all_transitions=False)
```

Loads and preprocesses reusable ExoMol, ExoAtom, ExoMolHR, HITRAN, or HITEMP
data. `auto` keeps small inputs in memory and converts large inputs to Parquet;
`parquet` forces persistent Parquet caching; `none` uses the original input
format. Use the result with `data=` in `px.conversion`,
`px.partition_functions`, `px.specific_heats`, `px.lifetimes`,
`px.cooling_functions`, `px.oscillator_strengths`, `px.stick_spectra`,
`px.cross_sections`, or `px.stick_spectra_cross_section`.

**Parameters**

| Parameter | Type | Required | Description |
|---|---|---|---|
| `inp_filepath` | `str` or `None` | No | Optional `.inp` file; keyword arguments override file values |
| `database` | `str` | Yes* | `ExoMol`, `ExoAtom`, `ExoMolHR`, `HITRAN`, or `HITEMP` |
| `molecule` | `str` | Conditional | Required for molecular databases |
| `atom` | `str` | Conditional | Required for ExoAtom |
| `isotopologue` | `str` | Conditional | Required except for ExoAtom |
| `dataset` | `str` | Yes* | Line-list dataset name |
| `species_id` | `int` or `None` | No | Species identifier; needed by conversions that write HITRAN IDs |
| `read_path` | `str` | Yes* | Database root containing the input line list |
| `save_path` | `str` | No | Output root retained for later calculations; default `./output/` |
| `logs_path` | `str` or `None` | No | Log path retained for later calculations; `None` disables file logging |
| `min_range` | `float` | No | Initial minimum wavenumber/wavelength; default `0` |
| `max_range` | `float` | No | Initial maximum wavenumber/wavelength; default `30000` |
| `wn_wl` | `str` | No | Range coordinate, `WN` or `WL`; default `WN` |
| `wn_wl_unit` | `str` | No | `cm-1`, `um`, or `nm`; default `cm-1` |
| `unc_filter` | `float` or `None` | No | Maximum uncertainty retained during preprocessing; default `None` |
| `cache` | `str` | No | `auto`, `parquet`, or `none`; default `auto` |
| `cache_dir` | `str` or `None` | No | Persistent cache directory; default is source-adjacent `.pyexocross_cache/` |
| `max_memory` | `int` | No | Auto-mode memory threshold in MB; default `512` |
| `refresh` | `bool` | No | Rebuild matching persistent caches; default `False` |
| `all_transitions` | `bool` | No | Eagerly load all transition files; default `False` |
| `device` or `run_mode` | `str` | No | `CPU` or `GPU`; default `CPU` |
| `gpu_backend` | `str` | No | GPU backend policy; default `AUTO` |
| `gpu_batch_lines` | `int` | No | GPU line-batch limit; default `8192` |
| `gpu_batch_grid` | `int` | No | GPU grid-batch limit; default `256` |
| `chunk_size` | `int` | No | Rows read per transition chunk; default `100000` |

`*` Not required when supplied by `inp_filepath`.

Later calculations may provide their own `min_range` and `max_range` as long
as that interval is covered by the data loaded here. A wider interval requires
another `px.load`. The `cutoff` value does not expand the range required from
`px.load`.
Range errors report the exact recommended `min_range` and
`max_range` for the next load.

A later `unc_filter` may keep the loaded value or use a smaller, stricter value;
PyExoCross then filters the loaded data further. A larger value, or changing a
numeric filter to `None`, requires another `px.load` and the error reports the
required filter value.

For ExoMol and ExoAtom, a later `qns_filter` may add constraints or narrow the
values retained by `px.load`. Empty value lists are wildcards and do not filter
the loaded states. Removing an active constraint or adding values excluded by
`px.load` requires loading the data again. HITRAN, HITEMP, and ExoMolHR retain
the complete line table and apply QN filters during each calculation.

For ExoMol/ExoAtom, calling `px.lifetimes`, `px.cooling_functions`, or
`px.oscillator_strengths` automatically expands range-loaded data to all
transitions and updates the same object for reuse. \
Set `all_transitions=True` only when this work should be done eagerly. \
Partition functions and specific heats reuse states and database metadata from
LoadedData but do not read or expand transitions.

**Example**

```python
data = px.load(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/MgH.log',
    min_range=0,
    max_range=10000,
    cache='auto',  
    max_memory=512,
    device='CPU',
)

# A narrower calculation range reuses the loaded data.
px.cross_section(
    data=data,
    min_range=1000,
    max_range=2000,
    temperatures=[296, 1000],
    pressures=[1.0],
)

# This first whole-list calculation automatically loads all transitions.
px.lifetimes(data=data)

# These calls reuse the automatically expanded object.
px.cooling_functions(data=data)
px.oscillator_strengths(data=data)
```

The original direct API remains available and still performs loading,
preprocessing, calculation, and saving in one call.

The direct API can also use the same cache without calling `px.load`:

```python
px.conversion(..., cache='auto')
px.lifetimes(..., cache='auto')
px.cooling_functions(..., cache='auto')
px.oscillator_strengths(..., cache='auto')
px.stick_spectra(..., cache='auto')
px.cross_section(..., cache='auto')
px.stick_spectra_cross_section(..., cache='auto')
```

With `cache_dir=None`, PyExoCross stores persistent files in
`<read_path>/<data_info>/.pyexocross_cache/`. Initial conversion takes extra
time; subsequent runs reuse the cache if the source path, size, and
modification time are unchanged. Set `refresh=True` to rebuild it.
`max_memory` is specified in MB. The default `max_memory=512` is converted
internally to `512 * 1024**2` bytes and only limits transition data retained by
`cache='auto'`; it is not a limit on total process memory.

The same Parquet files work in CPU and GPU modes. Current Parquet reading,
decompression, and preprocessing run on the CPU through PyArrow. GPU mode
accelerates supported calculation kernels after the required arrays are loaded.

For range-based spectra and cross sections, Parquet mode stores filtered files
under `<cache_dir>/ranges/` using names with specified ranges. \
A finite `cutoff` expands the cached interval only inside the transition files selected by the original
calculation range; it does not add adjacent source files. Cache validity and
source metadata are stored in `<cache_dir>/manifest.json`. When several valid
caches cover a request, PyExoCross selects the narrowest interval and prints
the number, intervals, and total size of range cache files. \
The complete states table is also stored in a parquet file. \
Its source path, size, and modification time
are tracked in the same manifest instead of being embedded in the filename.

HITRAN/HITEMP fixed-width files and ExoMolHR CSV files are normalized once and
stored as database-prefixed, range-named `*__00000-30000.linelist.parquet`
files. Their existing `v`
column is used for row-group range filtering, while uncertainty and intensity
filters remain calculation parameters. With `cache_dir=None`, these caches are
stored in `.pyexocross_cache` beside the source `.par` or `.csv` file.

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

***Note***

If you edit the same `.inp` file repeatedly in a long-lived session (e.g. Jupyter),
use `force_reload=True` to avoid reusing stale cached parameters.

**Example**

```python
px.run('/path/to/MgH_ExoMol.inp')

# Re-read updated .inp content in the same session
px.run('/path/to/MgH_ExoMol.inp', force_reload=True)
```

---

## `px.download`

```python
px.download(database, file_path=None, **kwargs)
```

Download database files from external servers. Supported databases are ExoMol, ExoMolHR, ExoAtom, and HITRAN.

**Parameters**

| Parameter | Type | Required | Description |
|---|---|---|---|
| `database` | `str` | Yes | Name of the database to download: `'ExoMol'`, `'ExoMolHR'`, `'ExoAtom'`, or `'HITRAN'`. |
| `file_path` | `str` | Yes | Path to save the downloaded files (e.g. `/path/to/database/`), can also use `save_path`. |
| `species_info` | `dict` | Yes | Molecule/isotopologue or atom configuration dictionary. |
| `download` | `bool` | Yes | If `True` (default), download the files. If `False`, only generate the URL file. |
| `wn_range` | `list` | Optional | Filter transitions by wavenumber range (cm⁻¹). |
| `T` | `int` | Optional | Temperature in Kelvin. |
| `threshold` | `float` | Optional | Minimum line intensity threshold in cm/molecule. |
| `dataset` | `str` | Optional | The dataset to use. Can be either `'NIST'` or `'Kurucz'`. |

### Species Configuration Options (`species_info`)

Inside `species_info`, you can configure the following optional settings:
- **`wn_range`** (list/tuple of `[wn_min, wn_max]`): Filter transitions by wavenumber range (cm⁻¹). Supported by ExoMol (downloads only segmented transition files fully within range), ExoMolHR, and HITRAN.
- **`T`** or **`temperature`** (int): Temperature in Kelvin. Required for ExoMolHR online query downloads.
- **`threshold`** or **`Smin`** (float): Minimum line intensity threshold in cm/molecule. Required for ExoMolHR online query downloads to filter out weaker lines.
- **`dataset`** (str): The dataset to use. Required for ExoAtom downloads. Can be either `'NIST'` (critically evaluated and recommended for accuracy) or `'Kurucz'` (recommended for completeness).

| Database | Required Parameters | Optional Parameters |
| :--- | :--- | :--- |
| `'ExoMol'` | `database`, `file_path`, `species_info`, `download` | `wn_range` |
| `'ExoMolHR'` | `database`, `file_path`, `species_info`, `download` | `T`, `threshold`, `wn_range` |
| `'ExoAtom'`| `database`, `file_path`, `species_info`, `download` | `dataset` |
| `'HITRAN'` | `database`, `file_path`, `species_info`, `download` | `wn_range` |

**Examples**

### ExoMol Download
```python
px.download(
    database='ExoMol',
    file_path='/path/to/ExoMol/',
    species_info={
        'MgH': {
            '24Mg-1H': {'wn_range': None},
            '25Mg-1H': {'wn_range': None},
        },
        'H2O': {
            '1H2-16O': {'wn_range': [41000, 41200]},
        },
    },
    download=True,
)
```

### ExoMolHR Download
```python
px.download(
    database='ExoMolHR',
    file_path='/path/to/ExoMolHR/',
    species_info={
        'MgH': {
            '24Mg-1H': None,
            '25Mg-1H': {'T': 1000, 'wn_range': [0, 500], 'threshold': 1e-30},
        },
        'AlH': {
            '27Al-1H': {'T': 500, 'wn_range': [0, 500], 'threshold': 1e-30}
        },
    },
    download=True,
)
```

### ExoAtom Download
```python
px.download(
    database='ExoAtom',
    file_path='/path/to/ExoAtom/',
    species_info={
        'He': {'dataset': 'NIST'},
        'He_p': {'3He_p': {'dataset': 'NIST'}},
        'Ar_p': {'dataset': 'Kurucz'},
    },
    download=True,
)
```

### HITRAN Download
```python
px.download(
    database='HITRAN',
    file_path='/path/to/HITRAN/',
    species_info={
        'NO': {
            '14N-16O': {'wn_range': [0, 100]},
            '15N-16O': {'wn_range': [100, 150]},
        },
        'H2O': {
            '1H2-16O': {'wn_range': [100, 110]},
        },
    },
    download=True,
)
```

### Generating URLs Only

If you only want to generate the list of URLs to be downloaded later, set `download=False`:

```python
px.download(
    database='ExoMol',
    file_path='/path/to/ExoMol/',
    species_info={'MgH': {'24Mg-1H': None}},
    download=False,
)
```

This will save the URLs to `/path/to/ExoMol/url/exomol__urls.txt` and print them to the console.

---

## `px.conversion`

```python
px.conversion(inp_filepath=None, **kwargs)
```

Convert between ExoMol/ExoMolHR/ExoAtom and HITRAN line-list formats.

**Parameters**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `inp_filepath` | `str` | `None` | Path to `.inp` file |
| `database` | `str` | `'ExoMol'` | `'ExoMol'`, `'ExoMolHR'`, `'ExoAtom'`, `'HITRAN'`, or `'HITEMP'` |
| `molecule` | `str` | `None` | Molecule name (e.g. `'MgH'`, `'NO'`) |
| `atom` | `str` | `None` | Atom name for ExoAtom (e.g. `'Ar'`, `'Li'`) |
| `isotopologue` | `str` | `None` | Isotopologue (e.g. `'24Mg-1H'`, `'14N-16O'`) |
| `dataset` | `str` | `None` | Dataset name (e.g. `'XAB'`, `'NIST'`) |
| `species_id` | `int` | `0` | Species identifier (e.g. `501`, `81`) |
| `read_path` | `str` | `'./'` | Path to input database |
| `save_path` | `str` | `'./output/'` | Path for output files |
| `logs_path` | `str` | `None` | Log file path |
| `conversion_format` | `str` | `None` | `'ExoMo'` = HITRAN/HITEMP -> ExoMol/ExoAtom, `'HITRAN'` = ExoMol/ExoMolHR/ExoAtom -> HITRAN/HITEMP |
| `conversion_min_freq` | `float` | `0` | Minimum frequency in cm⁻¹ |
| `conversion_max_freq` | `float` | `1e10` | Maximum frequency in cm⁻¹ |
| `conversion_unc` | `float` or `None` | `None` | Uncertainty filter (cm⁻¹); `None` to disable |
| `conversion_threshold` | `float` or `None` | `None` | Intensity threshold (cm/molecule); `None` to disable |
| `global_qn_label_list` | `list[str]` | `[]` | Global quantum number labels for HITRAN output |
| `global_qn_format_list` | `list[str]` | `[]` | Global QN format specifiers (total 15 chars in HITRAN2004) |
| `local_qn_label_list` | `list[str]` | `[]` | Local quantum number labels for HITRAN output |
| `local_qn_format_list` | `list[str]` | `[]` | Local QN format specifiers (total 15 chars in HITRAN2004) |
| `qnslabel_list` | `list[str]` | `[]` | Optional quantum number label override; ExoMol/ExoMolHR/ExoAtom labels are derived from definition metadata by default |
| `qnsformat_list` | `list[str]` | `[]` | Optional quantum number format override; ExoMol/ExoMolHR/ExoAtom formats are derived from definition metadata by default |
| `ncputrans` | `int` | `4` | CPU cores for processing transitions |
| `ncpufiles` | `int` | `1` | Files processed simultaneously |
| `chunk_size` | `int` | `100000` | Chunk size for transitions |
| `cache` | `str` | `'auto'` | `'auto'`, `'parquet'`, or `'none'` |
| `cache_dir` | `str` or `None` | `None` | Persistent cache directory |
| `max_memory` | `int` | `512` | Auto-mode memory threshold in MB |

**Legacy Aliases**

```python
px.convert_exomol_to_hitran(...)     # Sets conversion_format='HITRAN'
px.convert_exomolhr_to_hitran(...)   # Sets conversion_format='HITRAN'
px.convert_exoatom_to_hitran(...)    # Sets conversion_format='HITRAN'
px.convert_hitran_to_exomol(...)     # Sets conversion_format='ExoMol'
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
    conversion_format='HITRAN',
    global_qn_label_list=['ElecState', 'v', 'Omega'],
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

**Parameters**

Same as [`px.partition_functions`](#pxpartition_functions), plus optional
compute backend kwargs:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `device` or `run_mode` | `str` | `'CPU'` | `'CPU'` or `'GPU'` |
| `gpu_backend` | `str` | `'AUTO'` | `'AUTO'`, `'CUDA'`, `'PyTorch-CUDA'`, `'CuPy-CUDA'`, or `'MPS'` |
| `gpu_batch_lines` | `int` | `8192` | GPU line-batch size (memory control) |
| `gpu_batch_grid` | `int` | `256` | GPU grid-batch size (memory control) |

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
    device='CPU',
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
| `abundance` | `float` | `1.0` | Optional isotopic abundance multiplier |
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
| `qnslabel_list` | `list[str]` | `[]` | Optional QN label override; ExoMol/ExoAtom labels are derived from definition metadata by default |
| `qnsformat_list` | `list[str]` | `[]` | Optional QN format override; ExoMol/ExoAtom formats are derived from definition metadata by default |
| **Computing** | | | |
| `ncputrans` | `int` | `4` | CPU cores for transitions |
| `ncpufiles` | `int` | `1` | Files processed simultaneously |
| `chunk_size` | `int` | `100000` | Chunk size |
| `device` or `run_mode` | `str` | `'CPU'` | `'CPU'` or `'GPU'` |
| `gpu_backend` | `str` | `'AUTO'` | `'AUTO'`, `'CUDA'`, `'PyTorch-CUDA'`, `'CuPy-CUDA'`, or `'MPS'` |
| `gpu_batch_lines` | `int` | `8192` | GPU line-batch size (memory control) |
| `gpu_batch_grid` | `int` | `256` | GPU grid-batch size (memory control) |
| **Plotting** | | | |
| `plot` | `bool` | `False` | Whether to generate a plot |
| `plot_method` | `str` | `'log'` | `'log'` or `'linear'` |
| `plot_wn_wl` | `str` | `'WN'` | X-axis: `'WN'` or `'WL'` |
| `plot_unit` | `str` | `'cm-1'` | `'cm-1'`, `'um'`, or `'nm'` |
| `limit_yaxis` | `float` | `1e-30` | Lower limit for y-axis (cm/molecule) |

***Note***

`wn_wl` and `wn_wl_unit` define the calculation range, saved file name, and the
first column of the saved `.stick` file. \
Use `wn_wl='WN', wn_wl_unit='cm-1'`
to save wavenumber, or `wn_wl='WL', wn_wl_unit='nm'/'um'` to save wavelength. \
`plot_wn_wl` and `plot_unit` only control the plot x-axis.

### Quantum Number Filter (`qns_filter`)

The `qns_filter` parameter accepts a dictionary where keys are quantum number
labels and values are lists of accepted value patterns.  Use an empty list
`[]` to accept all values for that label.

```python
qns_filter={
    '+/-': [],                    # accept all parities
    'v': ['0,', '1,', '2,'],      # only v' = 0, 1, 2
    'ElecState': [],              # accept all electronic states
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
    abundance=1.0,
    device='CPU',
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

***Note***

The cross-section grid is built in the selected `wn_wl_unit`. \
For example,
with `wn_wl='WL', wn_wl_unit='nm'`, `min_range`, `max_range`, `bin_size`, the
saved file name, and the first column of the saved `.xsec` file all use nm. \
The plot axis can be changed independently with `plot_wn_wl` and `plot_unit`.
Line-profile widths and `cutoff` remain in cm⁻¹.


### Available Line Profiles

| Profile Name | Description |
|---|---|
| `'Doppler'` | Pure Doppler profile |
| `'Gaussian'` | Gaussian profile |
| `'Lorentzian'` | Lorentzian profile |
| `'SciPyVoigt'` | Voigt via SciPy `wofz` (recommended) |
| `'SciPyWofzVoigt'` | Voigt via SciPy `wofz` (alternative) |
| `'HumlicekVoigt'` | Humlicek algorithm for Voigt |
| `'PseudoVoigt'` | Generic pseudo-Voigt |
| `'PseudoThompsonVoigt'` | Thompson pseudo-Voigt approximation |
| `'PseudoKielkopfVoigt'` | Kielkopf pseudo-Voigt approximation |
| `'PseudoOliveroVoigt'` | Olivero pseudo-Voigt approximation |
| `'PseudoLiuLinVoigt'` | Liu-Lin pseudo-Voigt approximation |
| `'PseudoRoccoVoigt'` | Rocco pseudo-Voigt approximation |
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
| `<molecule>` or `<atom>` | <molecule> or <atom> broadening |

Support all broadening files format from ExoMol database.

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
    abs_emi='Ab',
    abundance=1,
    profile='SciPyVoigt',
    device='GPU',
    gpu_backend='AUTO',
    gpu_batch_lines=8192,
    gpu_batch_grid=256,
    broadeners=['Default'],
    ratios=[1.0],
    cutoff=25.0,
    plot=True,
    plot_method='log',
)
```

---

## `px.stick_spectra_cross_section`

```python
px.stick_spectra_cross_section(inp_filepath=None, **kwargs)
```

Calculate LTE or Non-LTE stick spectra and cross sections simultaneously.

**Parameters**

All parameters from [`px.stick_spectra`](#pxstick_spectra) and [`px.cross_sections`](#pxcross_sections).

The same `wn_wl`/`wn_wl_unit` selection is used for both saved outputs:
`.stick` and `.xsec` first columns are both wavenumber for `WN` and wavelength
for `WL`. Cross-section `bin_size` follows the same unit.

**Example**

```python
px.stick_spectra_cross_section(
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
    abs_emi='Ab',
    abundance=1,
    profile='SciPyVoigt',
    device='GPU',
    gpu_backend='AUTO',
    gpu_batch_lines=8192,
    gpu_batch_grid=256,
    plot=True,
)
```
