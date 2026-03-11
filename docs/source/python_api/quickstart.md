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

print(px.__version__)  # '1.0.0'
```

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
    logs_path='/path/to/output/log/exomol_xsec.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
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

## Minimal Example: Cross Sections from ExoMolHR

```python
import pyexocross as px

px.cross_sections(
    database='ExoMolHR',
    molecule='NO',
    isotopologue='14N-16O',
    read_path='/path/to/Databases/ExoMolHR/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomolhr_xsec.log',
    ncputrans=1,    
    ncpufiles=1,
    chunk_size=100000,
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
    logs_path='/path/to/output/log/exoatom_stick_nlte.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
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
    logs_path='/path/to/output/log/hitran_xsec.log',
    chunk_size=100000,
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
    broadeners=['Air', 'Self'],
    ratios=[0.7, 0.3],
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

### Logging

All API calls support a `logs_path` parameter.  When provided, screen output
is automatically duplicated to the log file.

```python
px.partition_functions(
    ...,
    logs_path='/path/to/output/log/my_run.log',
)
```

## Next Steps

- [**Functions Reference**](functions.md) -- Detailed documentation for every function
- [**Parameters Reference**](parameters.md) -- Complete list of all parameters with types and defaults
- [**Examples**](examples.md) -- Full working examples for ExoMol, ExoAtom, and HITRAN
- [**`.inp` File Mapping**](inp_mapping.md) -- Mapping between `.inp` keywords and Python kwargs
