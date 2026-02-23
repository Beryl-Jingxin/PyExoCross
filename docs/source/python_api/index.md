# Python API Reference

PyExoCross provides a full-featured Python API for programmatic access to all
functionalities.  Instead of writing `.inp` configuration files, you can call
functions directly from Python scripts, Jupyter notebooks, or any other Python
environment.

```python
import pyexocross as px
```

## Overview

| Function | Description |
|---|---|
| {func}`px.run <pyexocross.run>` | Run all enabled functions from an `.inp` file |
| {func}`px.conversion <pyexocross.conversion>` | Convert between ExoMol and HITRAN formats |
| {func}`px.partition_functions <pyexocross.partition_functions>` | Calculate partition functions $Q(T)$ |
| {func}`px.specific_heats <pyexocross.specific_heats>` | Calculate specific heats $C_p(T)$ |
| {func}`px.cooling_functions <pyexocross.cooling_functions>` | Calculate cooling functions $W(T)$ |
| {func}`px.lifetimes <pyexocross.lifetimes>` | Calculate radiative lifetimes $\tau$ |
| {func}`px.oscillator_strengths <pyexocross.oscillator_strengths>` | Calculate oscillator strengths $gf$ or $f$ |
| {func}`px.stick_spectra <pyexocross.stick_spectra>` | Calculate LTE / Non-LTE stick spectra |
| {func}`px.cross_sections <pyexocross.cross_sections>` | Calculate LTE / Non-LTE cross sections |

## Supported Databases

| Database | Keyword | Species Identifier |
|---|---|---|
| **ExoMol** | `database='ExoMol'` | `molecule`, `isotopologue`, `dataset` |
| **ExoAtom** | `database='ExoAtom'` | `atom`, `dataset` |
| **HITRAN** | `database='HITRAN'` | `molecule`, `isotopologue`, `dataset` |
| **HITEMP** | `database='HITEMP'` | `molecule`, `isotopologue`, `dataset` |

## Two Ways to Use PyExoCross

### 1. From an `.inp` file (traditional)

```python
px.run('/path/to/config.inp')
```

### 2. With keyword arguments (recommended for scripting)

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
    profile='SciPyVoigt',
)
```

You can also mix the two approaches -- pass an `.inp` file **and** override
individual parameters with keyword arguments:

```python
px.cross_sections(
    inp_filepath='/path/to/base_config.inp',
    temperatures=[3000],          # override .inp value
    pressures=[0.1, 1.0, 10.0],   # override .inp value
)
```

## Documentation Contents

```{toctree}
:maxdepth: 2

quickstart
functions
parameters
examples
inp_mapping
```
