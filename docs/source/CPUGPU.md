# CPU and GPU

This page explains how to run ***PyExoCross*** in CPU mode or GPU mode, and which functions support GPU acceleration.

## CPU mode (Default)

If you do not set GPU options, ***PyExoCross*** runs on CPU by default.

```bash
RunMode                                 CPU                       # CPU(default) or GPU
```

```python
import pyexocross as px

px.cross_sections(
    ...,
    run_mode='CPU',
    ...,
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

```bash
RunMode                                 GPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default), CUDA, PyTorch-CUDA, CuPy-CUDA, or MPS (used only when RunMode=GPU)
GPUBatchLines                           8192                      # GPU line-batch size (only used when RunMode=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when RunMode=GPU)
```

```python
import pyexocross as px

px.cross_sections(
    ...,
    run_mode='GPU',
    gpu_backend='AUTO',
    gpu_batch_lines=8192,
    gpu_batch_grid=256,
    ...,
)
```

## GPU backend options

| `gpu_backend` | Meaning | Selection / fallback order |
|---|---|---|
| `AUTO` | Recommended automatic mode | `PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU` |
| `CUDA` | CUDA policy mode | `PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU` |
| `PyTorch-CUDA` | Force PyTorch CUDA only | `PyTorch-CUDA -> CPU` |
| `CuPy-CUDA` | Force CuPy CUDA only | `CuPy-CUDA -> CPU` |
| `MPS` | Force Apple Metal backend | `MPS -> CPU` |

***Notes:***
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
