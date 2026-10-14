# CPU and GPU

This page explains how to run ***PyExoCross*** in CPU mode or GPU mode, and which functions support GPU acceleration.

## CPU mode (Default)

If you do not set GPU options, ***PyExoCross*** runs on CPU by default.

```bash
Device                                  CPU                       # CPU(default) or GPU
```

```python
import pyexocross as px

px.cross_sections(
    ...,
    device='CPU',
    ...,
)
```

Use CPU mode when:
- no GPU is available
- you want maximum portability
- you are running small jobs where GPU startup overhead is unnecessary

## GPU mode (CUDA or MPS)

Enable GPU mode with:
- `device='GPU'`
- optional `gpu_backend` selection

```bash
Device                                  GPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default), CUDA, PyTorch-CUDA, CuPy-CUDA, or MPS (used only when Device=GPU)
GPUBatchLines                           8192                      # GPU line-batch size (only used when Device=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when Device=GPU)
```

```python
import pyexocross as px

px.cross_sections(
    ...,
    device='GPU',
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

- If the requested backend is not available, ***PyExoCross*** falls back to CPU.

## GPU batch sizes

- `gpu_batch_lines` is the maximum number of spectral lines processed in one GPU batch.
- `gpu_batch_grid` is the maximum number of spectral grid points processed in one GPU batch.

The main temporary arrays have an approximate shape of
`(gpu_batch_grid, gpu_batch_lines)`, so GPU memory use generally increases
with the product of the two values. Larger batches may reduce launch overhead
but require more GPU memory. If an out-of-memory error occurs, reduce one value
at a time, for example by a factor of two.

Both parameters may be any positive integer. They are data-batch dimensions,
not CUDA thread-block dimensions, and therefore do not need to be powers of two
or multiples of 32. The final incomplete batch is handled automatically.

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
