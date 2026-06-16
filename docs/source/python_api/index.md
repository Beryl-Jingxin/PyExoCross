# <font color=Orange>**Python API Reference**</font>

PyExoCross provides a full-featured Python API for programmatic access to all
functionalities.  Instead of writing `.inp` configuration files, you can call
functions directly from Python scripts, Jupyter notebooks, or any other Python
environment.

```bash
pip install pyexocross
```

```python
import pyexocross as px
```

## <font color=Orange>**PyExoCross Documentation Homepage**</font>

Find more details of input parameters from [***PyExoCross Documentation Homepage***](../index.rst).

## Overview

| Function | Description |
|---|---|
| {func}`px.run <pyexocross.run>` | Run all enabled functions from an `.inp` file |
| {func}`px.conversion <pyexocross.conversion>` | Convert between ExoMol and HITRAN data formats |
| {func}`px.partition_functions <pyexocross.partition_functions>` | Calculate partition functions $Q(T)$ |
| {func}`px.specific_heats <pyexocross.specific_heats>` | Calculate specific heats $C_p(T)$ |
| {func}`px.cooling_functions <pyexocross.cooling_functions>` | Calculate cooling functions $W(T)$ |
| {func}`px.lifetimes <pyexocross.lifetimes>` | Calculate radiative lifetimes $\tau$ |
| {func}`px.oscillator_strengths <pyexocross.oscillator_strengths>` | Calculate oscillator strengths $gf$ or $f$ |
| {func}`px.stick_spectra <pyexocross.stick_spectra>` | Calculate LTE / Non-LTE stick spectra |
| {func}`px.cross_sections <pyexocross.cross_sections>` | Calculate LTE / Non-LTE cross sections |
| {func}`px.stick_spectra_cross_section <pyexocross.stick_spectra_cross_section>` | Calculate LTE / Non-LTE stick spectra and cross sections simultaneously |
| {func}`px.download <pyexocross.download>` | Download database files (ExoMol, ExoAtom, ExoMolHR, HITRAN) |

## Supported Databases

| Database | Keyword | Species Identifier |
|---|---|---|
| **ExoMol**   | `database='ExoMol'`   | `molecule`, `isotopologue`, `dataset` |
| **ExoMolHR** | `database='ExoMolHR'` | `molecule`, `isotopologue` |
| **ExoAtom**  | `database='ExoAtom'`  | `atom`, `dataset` |
| **HITRAN**   | `database='HITRAN'`   | `molecule`, `isotopologue`, `dataset` |
| **HITEMP**   | `database='HITEMP'`   | `molecule`, `isotopologue`, `dataset` |

## Two Ways to Use PyExoCross Python Package

### 1. From an `.inp` file (traditional)

```python
px.run('/path/to/config.inp')

# Force re-parse when the same .inp file was edited in this session
px.run('/path/to/config.inp', force_reload=True)
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

## CPU / GPU Compute Mode

PyExoCross uses CPU mode by default.  You can switch to GPU mode with
`run_mode='GPU'` and optionally select a backend with `gpu_backend`.

```python
# Default CPU mode
px.cross_sections(..., run_mode='CPU')

# Auto-select GPU backend (recommended)
# Priority: PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback
px.cross_sections(..., run_mode='GPU', gpu_backend='AUTO')

# CUDA policy
# Priority: PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback
px.cross_sections(..., run_mode='GPU', gpu_backend='CUDA')

# Force PyTorch CUDA only
px.cross_sections(..., run_mode='GPU', gpu_backend='PyTorch-CUDA')

# Force CuPy CUDA only
px.cross_sections(..., run_mode='GPU', gpu_backend='CuPy-CUDA')

# Force Apple Metal (MPS)
px.cross_sections(..., run_mode='GPU', gpu_backend='MPS')
```

GPU memory-control knobs:

```python
px.cross_sections(
    ...,
    run_mode='GPU',
    gpu_backend='AUTO',   # AUTO (recommended), CUDA, PyTorch-CUDA, CuPy-CUDA, or MPS
    gpu_batch_lines=8192,
    gpu_batch_grid=256,
)
```

***Notes:***
- If GPU backend is unavailable, execution falls back to CPU formulas.
- `AUTO` order: `PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU`.
- `CUDA` order: `PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU`.
- On Apple Silicon, MPS uses float32 kernels so tiny numeric differences vs CPU/CUDA are expected.
- GPU acceleration applies to:
  `cooling_functions`, `stick_spectra`, `cross_sections`,
  `stick_spectra_cross_section`.

You can also mix the two approaches -- pass an `.inp` file **and** override
individual parameters with keyword arguments:

```python
px.cross_sections(
    inp_filepath='/path/to/base_config.inp',
    temperatures=[3000],          # override .inp value
    pressures=[0.1, 1.0, 10.0],   # override .inp value
)
```

## PyExoCross Python Package Documentation Contents

```{toctree}
:maxdepth: 2

quickstart
functions
parameters
examples
inp_mapping
```
