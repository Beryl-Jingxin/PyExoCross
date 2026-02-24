# Quick Start

This guide walks you through the basics of using PyExoCross as a Python
library in under 5 minutes.

## Installation

```bash
# Clone the repository
git clone https://github.com/ExoMol/PyExoCross.git
cd PyExoCross

# Install dependencies
pip install -r requirements.txt

# Install PyExoCross in development mode
pip install -e .
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
    species_id=501,
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/run.log',
    temperatures=[1000],
    pressures=[1.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    bin_size=0.1,
    profile='SciPyVoigt',
    broadeners=['Default'],
    ratios=[1.0],
)
```

## Minimal Example: Cross Sections from ExoAtom

```python
import pyexocross as px

px.cross_sections(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/run.log',
    temperatures=[1000],
    pressures=[1.0],
    profile='SciPyVoigt',
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
    logs_path='/path/to/output/log/run.log',
    temperatures=[1000],
    pressures=[1.0],
    profile='SciPyVoigt',
)
```

## Using an `.inp` File

If you already have a `.inp` configuration file, you can run it directly:

```python
px.run('/path/to/MgH_ExoMol.inp')
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
| ExoMol | MgH (24Mg-1H) | `501` | HITRAN molecule ID 50, isotopologue ID 1 |
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
