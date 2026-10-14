# Quick Start

This guide walks you through the basics of using PyExoCross as a Python
library in under 5 minutes.

## Installation

```bash
# Install PyExoCross python package
pip install pyexocross
```

## Import

```python
import pyexocross as px

print(px.__version__)
```

## Choose CPU or GPU

PyExoCross defaults to CPU mode.

```python
# Explicit CPU mode (default)
run_mode='CPU'

# GPU mode with automatic backend selection (recommended)
run_mode='GPU'
gpu_backend='AUTO'   # PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback
gpu_batch_lines=8192
gpu_batch_grid=256

# CUDA mode (PyTorch-CUDA first, then CuPy-CUDA, then MPS fallback)
run_mode='GPU'
gpu_backend='CUDA'

# Force PyTorch CUDA only
run_mode='GPU'
gpu_backend='PyTorch-CUDA'

# Force CuPy CUDA only
run_mode='GPU'
gpu_backend='CuPy-CUDA'

# Force MPS (Apple Metal)
run_mode='GPU'
gpu_backend='MPS'
```

`gpu_batch_lines` and `gpu_batch_grid` set the maximum numbers of spectral
lines and spectral grid points processed in one GPU batch, respectively. GPU
memory use generally increases with their product. Both may be any positive
integer; powers of two and multiples of 32 are not required.

If no selected backend is available, PyExoCross falls back to CPU formulas.

GPU acceleration is available for `cooling_functions`, `stick_spectra`,
`cross_sections`, and `stick_spectra_cross_section`.  

## Minimal Example: Cross Sections from ExoMol

```python
import pyexocross as px

px.cross_sections(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_xsec.log',  # Use log='none' to disable file logging
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    run_mode='CPU',
    temperatures=[1000],
    pressures=[1.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    bin_size=0.1,
    profile='SciPyVoigt',
    predissociation=False,
    broadeners=['Default'],
    ratios=[1.0],
    cutoff=25.0,
    abs_emi='Ab',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN',
    plot_unit='cm-1',
    limit_yaxis=1e-30
)
```

## Keep the result in Python

The default `output='files'` saves files and returns `None`. For interactive
work, request an in-memory result:

```python
result = px.cross_sections(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    temperatures=[300],
    pressures=[1.0],
    min_range=0,
    max_range=100,
    bin_size=0.1,
    profile='Gaussian',
    output='memory',
    log='none',
)

xsec = result.select(T=300, P=1.0)
grid = result.coords['wavenumber']
```

Use `output='both'` to return the result and save files. See
[**Interactive results**](interactive_results.md) for combined calculations,
Non-LTE selection, units, parameters, timing, and help.

## Minimal Example: Cross Sections from ExoMolHR

```python
import pyexocross as px

px.cross_sections(
    database='ExoMolHR',
    molecule='NO',
    isotopologue='14N-16O',
    read_path='/path/to/Databases/ExoMolHR/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomolhr_xsec.log',  # Use log='none' to disable file logging
    ncputrans=1,    
    ncpufiles=1,
    chunk_size=100000,
    run_mode='CPU',
    temperatures=[1000],
    pressures=[1.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=24,
    max_range=53452,
    bin_size=0.1,
    profile='SciPyVoigt',
    predissociation=False,
    broadeners=['Default'],
    ratios=[1.0],
    alpha_hwhm=3.0,
    gamma_hwhm=None,
    cutoff=25.0,
    abs_emi='Ab',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN',
    plot_unit='cm-1',
    limit_yaxis=1e-30
)
```

## Minimal Example: Cross Sections from ExoAtom

```python
import pyexocross as px

px.cross_sections(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exoatom_stick_nlte.log',  # Use log='none' to disable file logging
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    run_mode='CPU',
    temperatures=[1000],
    pressures=[1.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=115400,
    bin_size=1,
    profile='SciPyVoigt',
    predissociation=False,
    broadeners=['Default'],
    ratios=[1.0],
    cutoff=25.0,
    abs_emi='Ab',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN',
    plot_unit='cm-1',
    limit_yaxis=1e-30
)
```

## Minimal Example: Cross Sections from HITRAN

```python
import pyexocross as px

px.cross_sections(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran_xsec.log',  # Use log='none' to disable file logging
    chunk_size=100000,
    run_mode='CPU',
    qnslabel_list=['J', 'X', 'Omega', 'v1', 'Sym', 'F'],
    qnsformat_list=['%5.1f', '%2s', '%3s', '%2d', '%1s', '%5s'],
    temperatures=[1000],
    pressures=[1.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    bin_size=0.1,
    profile='SciPyVoigt',
    broadeners=['Air'],
    ratios=[1.0],
    cutoff=25.0,
    abs_emi='Ab',
    plot=True,
    plot_method='log',
    plot_unit='cm-1',
    limit_yaxis=1e-30
)
```

## Using an `.inp` File

If you already have a `.inp` configuration file, you can run it directly:

```python
px.run('/path/to/MgH_ExoMol.inp')

# Use force_reload=True if you edited the same .inp file in this session
px.run('/path/to/MgH_ExoMol.inp', force_reload=True)
```

Or call a specific function from a `.inp` file (ignoring the on/off switches
in the file):

```python
px.cross_sections(inp_filepath='/path/to/MgH_ExoMol.inp')
```

## Key Concepts

### `read_path` -- Where Your Data Lives

- **ExoMol**: Path to the root ExoMol database directory. \
  Files are located as `{read_path}/{molecule}/{isotopologue}/{dataset}/`. \
  If a directory is given, the code looks for definition, line list, and partition function files in fowllowing format filepath \
  `{read_path}/{molecule}/{isotopologue}/{dataset}/{isotopologue}__{dataset}.def.json`. 
- **ExoMolHR**: Path to the root ExoMolHR database directory. \
  Files are located as `{read_path}/{molecule}/{isotopologue}/`. \
  If a directory is given, the code looks for line list and partition function files in fowllowing format filepath \
  `{read_path}/{molecule}/{isotopologue}/{isotopologue}__{dataset}.pf`.   
- **ExoAtom**: Path to the root ExoAtom database directory. \
  Files are located as `{read_path}/{atom}/{dataset}/`. \
  If a directory is given, the code looks for definition, line list, and partition function files in fowllowing format filepath \
  `{read_path}/{atom}/{dataset}/{atom}__{dataset}.adef.json`. 
- **HITRAN / HITEMP**: Path to the root HITRAN/HITEMP database directory. \
  Files are located as `{read_path}/{molecule}/{isotopologue}/`. \
  If a directory is given, the code looks for line list and partition function files in fowllowing format filepath \
  `{read_path}/{molecule}/{isotopologue}/{molecule}__{isotopologue}.par`.

### `save_path` -- Where Results Go

All output files (spectra, partition functions, plots, etc.) are saved under
`save_path`, organized into sub-directories by function type.

### `species_id` -- Species Identifier

A numeric identifier for the molecular or atomic species:

| Database | Example Species | `species_id` | Meaning |
|---|---|---|---|
| ExoMol | MgH (24Mg-1H) | `502` | HITRAN molecule ID 50, isotopologue ID 2 |
| ExoMolHR | NO (14N-16O) | `81` | HITRAN molecule ID 8, isotopologue ID 1 |
| ExoAtom | Ar | `601` | HITRAN atom ID 60, isotope ID 1 |
| HITRAN | NO (14N-16O) | `81` | HITRAN molecule ID 8, isotopologue ID 1 |

The first digits (`species_id // 10`) are the *species main ID*, and the last
digit (`species_id % 10`) is the *species sub ID*.  For ExoMol/HITRAN these
correspond to the HITRAN molecule and isotopologue numbers respectively.

### Logging and terminal output

All API calls support `logs_path`, `log`, and `verbose`. A valid `logs_path`
with `log='auto'` writes a log. Use `log='none'` to disable file logging.
Use `verbose=False` to hide normal terminal output and progress bars while
retaining any enabled log file.

```python
px.partition_functions(
    ...,
    logs_path='/path/to/output/log/my_run.log',  # Use log='none' to disable file logging
    log='auto',
    verbose=False,
)
```

For `.inp` workflows, use:

```text
LogFilePath                             None
Verbose                                 True
```

This shows normal terminal output without creating a log file.

## Next Steps

- [**Interactive results**](interactive_results.md) -- Use calculation results in Jupyter, Spyder, or Python
- [**Functions Reference**](functions.md) -- Detailed documentation for every function
- [**Parameters Reference**](parameters.md) -- Complete list of all parameters with types and defaults
- [**Examples**](examples.md) -- Full working examples for ExoMol, ExoAtom, and HITRAN
- [**`.inp` File Mapping**](inp_mapping.md) -- Mapping between `.inp` keywords and Python kwargs
