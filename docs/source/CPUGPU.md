# CPU / GPU Usage

This page explains how to run ***PyExoCross*** in CPU mode or GPU mode, and which functions support GPU acceleration.

## CPU mode (Default)

If you do not set GPU options, ***PyExoCross*** runs on CPU by default.

```python
import ***PyExoCross*** as px

px.cross_sections(
    ...,
    run_mode='CPU',
)
```

Use CPU mode when:
- no GPU is available
- you want maximum portability
- you are running small jobs where GPU startup overhead is unnecessary

## GPU mode (CUDA or MPS)

Enable GPU mode with:
- `run_mode='GPU'`
- optional `gpu_backend` selection

```python
import pyexocross as px

px.cross_sections(
    ...,
    run_mode='GPU',
    gpu_backend='AUTO',
    gpu_batch_lines=8192,
    gpu_batch_grid=256,
)
```

## GPU backend options

| `gpu_backend` | Meaning | Selection / fallback order |
|---|---|---|
| `AUTO` | Recommended automatic mode | `CuPy-CUDA -> PyTorch-CUDA -> MPS -> CPU` |
| `CUDA` | CUDA policy mode | `CuPy-CUDA -> PyTorch-CUDA -> MPS -> CPU` |
| `CuPy-CUDA` | Force CuPy CUDA only | `CuPy-CUDA -> CPU` |
| `PyTorch-CUDA` | Force PyTorch CUDA only | `PyTorch-CUDA -> CPU` |
| `MPS` | Force Apple Metal backend | `MPS -> CPU` |

Notes:
- `gpu_batch_lines` and `gpu_batch_grid` are memory-control knobs for GPU calculations.
- If the requested backend is not available, ***PyExoCross*** falls back to CPU.

## Functions support GPU

Only the following three functions support GPU acceleration:

| Function               | GPU Support |
|------------------------|-------------|
| Convert Data Format    | ❌          |
| Partition Functions    | ❌          |
| Specific Heats         | ❌          |
| Cooling Functions      | ✅          |
| Lifetimes              | ❌          |
| Oscillator Strengths   | ❌          |
| Stick Spectra          | ✅          |
| Cross Sections         | ✅          |

All other functions run on CPU.

