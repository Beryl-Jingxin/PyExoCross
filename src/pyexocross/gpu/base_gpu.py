"""GPU backend runtime selection for CUDA / MPS.

This module centralizes compute backend selection for PyExoCross:
- CPU
- GPU(CUDA / MPS)

CUDA can run with Torch or CuPy backend.
MPS uses Torch backend.
"""

from __future__ import annotations

from typing import Any, Optional, Tuple

import os
import shutil
import numpy as np
from tabulate import tabulate
import subprocess

from pyexocross.base.log import print_cpu_device_info

_DEFAULT_GPU_BATCH_LINES = 8192
_DEFAULT_GPU_BATCH_GRID = 256

_RUNTIME = {
    'requested_mode': 'CPU',
    'requested_gpu_backend': 'AUTO',
    'active_mode': 'CPU',
    'active_gpu_backend': 'CPU',
    'provider': None,  # 'torch' | 'cupy' | None
    'torch': None,
    'torch_checked': False,
    'cupy': None,
    'cupy_checked': False,
    'torch_device': None,
    'reason': 'GPU mode not requested',
    'gpu_batch_lines': _DEFAULT_GPU_BATCH_LINES,
    'gpu_batch_grid': _DEFAULT_GPU_BATCH_GRID,
}


def _gpu_debug_enabled() -> bool:
    return str(os.environ.get('PYEXOCROSS_GPU_DEBUG', '')).strip().lower() in ('1', 'true', 'yes', 'on')


def _drop_path_entries(path_value: str, blocked_tokens: tuple[str, ...]) -> str:
    if not path_value:
        return path_value
    kept = []
    for p in path_value.split(':'):
        p_l = p.lower()
        if any(tok in p_l for tok in blocked_tokens):
            continue
        kept.append(p)
    return ':'.join(kept)


def _detect_cuda_root_from_nvcc() -> Optional[str]:
    nvcc = shutil.which('nvcc')
    if not nvcc:
        return None
    real_nvcc = os.path.realpath(nvcc)
    cuda_root = os.path.dirname(os.path.dirname(real_nvcc))
    if os.path.isdir(cuda_root):
        return cuda_root
    return None


def _looks_like_rocm_path(value: Optional[str]) -> bool:
    if not value:
        return False
    v = value.lower()
    return ('/rocm' in v) or ('/hip' in v)


def _sanitize_cupy_cuda_env() -> None:
    """
    Sanitize process env for CuPy CUDA JIT compilation.

    In mixed ROCm/CUDA clusters, HIP include/library paths can leak into
    nvcc/nvrtc compile steps and break CuPy kernels (e.g. hip_fp16.h errors).
    """
    blocked = ('/rocm/', '/hip/')
    for key in ('PATH', 'LD_LIBRARY_PATH', 'CPATH', 'C_INCLUDE_PATH', 'CPLUS_INCLUDE_PATH', 'LIBRARY_PATH'):
        os.environ[key] = _drop_path_entries(os.environ.get(key, ''), blocked)

    # Remove explicit ROCm markers that may alter backend toolchain resolution.
    for key in list(os.environ.keys()):
        ku = key.upper()
        if ku.startswith('HIP_') or ku.startswith('ROCM_'):
            os.environ.pop(key, None)

    # Hint CUDA root explicitly based on current toolchain discovery.
    cuda_root = _detect_cuda_root_from_nvcc()
    if cuda_root:
        os.environ['CUDA_HOME'] = cuda_root
        os.environ['CUDA_PATH'] = cuda_root
        os.environ['CUPY_CUDA_PATH'] = cuda_root
        cxx = os.environ.get('CUDACXX')
        detected_nvcc = os.path.join(cuda_root, 'bin', 'nvcc')
        if os.path.isfile(detected_nvcc) and ((not cxx) or _looks_like_rocm_path(cxx)):
            os.environ['CUDACXX'] = detected_nvcc
            os.environ['NVCC'] = detected_nvcc


def _patch_cupy_cuda_path(cp: Any) -> None:
    """
    Override CuPy's cached CUDA toolchain paths when they point to ROCm.
    """
    cuda_root = os.environ.get('CUPY_CUDA_PATH') or os.environ.get('CUDA_PATH') or os.environ.get('CUDA_HOME')
    if not cuda_root:
        return
    if _looks_like_rocm_path(cuda_root):
        return
    try:
        import cupy._environment as cp_env  # type: ignore
    except Exception:
        return
    try:
        cp_env._cuda_path = cuda_root
        nvcc = os.path.join(cuda_root, 'bin', 'nvcc')
        if os.path.isfile(nvcc):
            cp_env._nvcc_path = nvcc
        # Clear ROCm/HIP toolchain hints for CUDA runs.
        if hasattr(cp_env, '_rocm_path'):
            cp_env._rocm_path = None
        if hasattr(cp_env, '_hipcc_path'):
            cp_env._hipcc_path = None
    except Exception:
        pass


def normalize_run_mode(value: Any) -> str:
    """Normalize run mode to ``CPU`` or ``GPU``."""
    if value is None:
        return 'CPU'
    mode = str(value).strip().upper()
    if mode in ('CPU', 'GPU'):
        return mode
    raise ValueError("run_mode must be 'CPU' or 'GPU'.")


def normalize_gpu_backend(value: Any) -> str:
    """Normalize GPU backend selection value.

    Accepted values:
    - AUTO
    - CUDA (legacy auto-select: prefer PyTorch then CuPy)
    - PYTORCH-CUDA
    - CUPY-CUDA
    - MPS
    """
    if value is None:
        return 'AUTO'
    backend = str(value).strip().upper()
    aliases = {
        'TORCH': 'PYTORCH-CUDA',
        'PYTORCH': 'PYTORCH-CUDA',
        'TORCHCUDA': 'PYTORCH-CUDA',
        'TORCH-CUDA': 'PYTORCH-CUDA',
        'PYTORCHCUDA': 'PYTORCH-CUDA',
        'CUPY': 'CUPY-CUDA',
        'CUPYCUDA': 'CUPY-CUDA',
    }
    backend = aliases.get(backend, backend)
    if backend in ('AUTO', 'CUDA', 'PYTORCH-CUDA', 'CUPY-CUDA', 'MPS'):
        return backend
    raise ValueError(
        "gpu_backend must be one of: 'AUTO', 'CUDA', 'PYTORCH-CUDA', 'CUPY-CUDA', or 'MPS'."
    )


def _normalize_positive_int(value: Any, name: str, default: int) -> int:
    if value in (None, 'None'):
        return default
    try:
        out = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a positive integer.") from exc
    if out <= 0:
        raise ValueError(f"{name} must be a positive integer.")
    return out


def _try_import_cupy() -> Optional[Any]:
    if _RUNTIME['cupy_checked']:
        return _RUNTIME['cupy']
    _RUNTIME['cupy_checked'] = True
    _sanitize_cupy_cuda_env()
    try:
        import cupy as cp  # type: ignore
    except Exception:
        _RUNTIME['cupy'] = None
        return None
    _patch_cupy_cuda_path(cp)
    _RUNTIME['cupy'] = cp
    return cp


def _try_import_torch() -> Optional[Any]:
    if _RUNTIME['torch_checked']:
        return _RUNTIME['torch']
    _RUNTIME['torch_checked'] = True
    try:
        import torch  # type: ignore
    except Exception:
        _RUNTIME['torch'] = None
        return None
    _RUNTIME['torch'] = torch
    return torch


def _select_cuda_provider() -> tuple[Optional[str], str, Optional[str], Optional[str]]:
    """Return (provider, reason, torch_device, active_backend)."""
    torch = _try_import_torch()
    if torch is not None and bool(torch.cuda.is_available()):
        return ('torch', 'CUDA enabled with Torch', 'cuda', 'CUDA')

    cp = _try_import_cupy()
    if cp is not None:
        try:
            n_dev = int(cp.cuda.runtime.getDeviceCount())
        except Exception:
            n_dev = 0
        if n_dev > 0:
            return ('cupy', f'CUDA enabled with CuPy ({n_dev} device(s))', None, 'CUDA')

    return (None, 'CUDA requested but no CUDA-capable runtime was detected', None, None)


def _select_torch_cuda_provider() -> tuple[Optional[str], str, Optional[str], Optional[str]]:
    """Return (provider, reason, torch_device, active_backend) for PyTorch-CUDA."""
    torch = _try_import_torch()
    if torch is not None and bool(torch.cuda.is_available()):
        return ('torch', 'CUDA enabled with Torch', 'cuda', 'CUDA')
    return (None, 'PyTorch-CUDA requested but Torch CUDA runtime is unavailable', None, None)


def _select_cupy_provider() -> tuple[Optional[str], str, Optional[str], Optional[str]]:
    """Return (provider, reason, torch_device, active_backend) for CuPy-CUDA."""
    cp = _try_import_cupy()
    if cp is not None:
        try:
            n_dev = int(cp.cuda.runtime.getDeviceCount())
        except Exception:
            n_dev = 0
        if n_dev > 0:
            return ('cupy', f'CUDA enabled with CuPy ({n_dev} device(s))', None, 'CUDA')
    return (None, 'CuPy-CUDA requested but CuPy CUDA runtime is unavailable', None, None)


def _select_mps_provider() -> tuple[Optional[str], str, Optional[str], Optional[str]]:
    """Return (provider, reason, torch_device, active_backend)."""
    torch = _try_import_torch()
    if torch is None:
        return (None, 'MPS requested but Torch is not available', None, None)
    try:
        mps_backend = getattr(torch.backends, 'mps', None)
        mps_ok = bool(mps_backend is not None and mps_backend.is_available())
    except Exception:
        mps_ok = False
    if mps_ok:
        return ('torch', 'MPS enabled with Torch', 'mps', 'MPS')
    return (None, 'MPS requested but no MPS-capable runtime was detected', None, None)


def configure_runtime(
    run_mode: Any = 'CPU',
    gpu_backend: Any = 'AUTO',
    gpu_batch_lines: Any = _DEFAULT_GPU_BATCH_LINES,
    gpu_batch_grid: Any = _DEFAULT_GPU_BATCH_GRID,
    verbose: bool = True,
) -> dict:
    """Configure global runtime backend."""
    requested_mode = normalize_run_mode(run_mode)
    requested_gpu_backend = normalize_gpu_backend(gpu_backend)

    batch_lines = _normalize_positive_int(gpu_batch_lines, 'gpu_batch_lines', _DEFAULT_GPU_BATCH_LINES)
    batch_grid = _normalize_positive_int(gpu_batch_grid, 'gpu_batch_grid', _DEFAULT_GPU_BATCH_GRID)

    _RUNTIME['requested_mode'] = requested_mode
    _RUNTIME['requested_gpu_backend'] = requested_gpu_backend
    _RUNTIME['gpu_batch_lines'] = batch_lines
    _RUNTIME['gpu_batch_grid'] = batch_grid

    # reset active
    _RUNTIME['active_mode'] = 'CPU'
    _RUNTIME['active_gpu_backend'] = 'CPU'
    _RUNTIME['provider'] = None
    _RUNTIME['torch_device'] = None

    if requested_mode == 'CPU':
        _RUNTIME['reason'] = 'GPU mode not requested'
    else:
        provider = None
        reason = ''
        torch_device = None
        active_backend = None

        if requested_gpu_backend in ('AUTO', 'CUDA'):
            provider, reason, torch_device, active_backend = _select_cuda_provider()
        elif requested_gpu_backend == 'PYTORCH-CUDA':
            provider, reason, torch_device, active_backend = _select_torch_cuda_provider()
        elif requested_gpu_backend == 'CUPY-CUDA':
            provider, reason, torch_device, active_backend = _select_cupy_provider()

        # Keep a CUDA->MPS fallback to support environments where input/API
        # still requests CUDA by default but only local MPS is available.
        if provider is None and requested_gpu_backend in ('AUTO', 'CUDA', 'MPS'):
            mps_provider, mps_reason, mps_device, mps_backend = _select_mps_provider()
            if mps_provider is not None:
                provider = mps_provider
                torch_device = mps_device
                active_backend = mps_backend
                if requested_gpu_backend == 'CUDA':
                    reason = 'CUDA requested but unavailable; fallback to MPS with Torch'
                else:
                    reason = mps_reason
            elif requested_gpu_backend == 'AUTO' and reason:
                reason = f'{reason} {mps_reason}'
            elif not reason:
                reason = mps_reason

        if provider is not None:
            _RUNTIME['active_mode'] = 'GPU'
            _RUNTIME['active_gpu_backend'] = active_backend if active_backend is not None else 'CPU'
            _RUNTIME['provider'] = provider
            _RUNTIME['torch_device'] = torch_device
            _RUNTIME['reason'] = reason
        else:
            _RUNTIME['reason'] = reason if reason else 'GPU backend unavailable; fallback to CPU'

    if verbose:
        if _RUNTIME["requested_mode"] == "GPU":
            print("\n=== NVIDIA-SMI ===")
            try:
                result = subprocess.run(["nvidia-smi"], capture_output=True, text=True, check=True)
                print(result.stdout)
            except Exception as exc:
                print("nvidia-smi unavailable:", exc)

            print("\n=== PyTorch CUDA ===")
            try:
                torch = _try_import_torch()
                if torch is None:
                    raise RuntimeError("PyTorch not installed")
                print("PyTorch:", getattr(torch, "__version__", "unknown"))
                print("CUDA available:", bool(torch.cuda.is_available()))
                print("CUDA build version:", getattr(torch.version, "cuda", None))
                if torch.cuda.is_available():
                    print("GPU count:", torch.cuda.device_count())
                    for i in range(torch.cuda.device_count()):
                        props = torch.cuda.get_device_properties(i)
                        print(f"\nGPU {i}: {torch.cuda.get_device_name(i)}")
                        print("  Total memory (GB):", props.total_memory / 1024**3)
                        print("  Compute capability:", f"{props.major}.{props.minor}")
                        print("  SM count:", props.multi_processor_count)
            except Exception as exc:
                print("PyTorch check failed:", exc)

            print("\n=== CuPy CUDA ===")
            try:
                cp = _try_import_cupy()
                if cp is None:
                    raise RuntimeError("CuPy not installed")
                n_dev = int(cp.cuda.runtime.getDeviceCount())
                print("CuPy:", getattr(cp, "__version__", "unknown"))
                print("CUDA devices:", n_dev)
                if _gpu_debug_enabled():
                    try:
                        print("CuPy CUDA path:", cp.cuda.get_cuda_path())
                    except Exception:
                        pass
                    print("CUDA_HOME:", os.environ.get("CUDA_HOME", ""))
                    print("CUDA_PATH:", os.environ.get("CUDA_PATH", ""))
                    print("CUDACXX:", os.environ.get("CUDACXX", ""))
            except Exception as exc:
                print("CuPy check failed:", exc)

        print_cpu_device_info(prefix="Host CPU")

        if _RUNTIME["requested_mode"] == "CPU":
            rows = [
                ["Requested mode", _RUNTIME["requested_mode"]],
                ["Active mode", _RUNTIME["active_mode"]],
                ["Provider", _RUNTIME["provider"]],
                ["Reason", _RUNTIME["reason"]],
            ]
        else:
            # Keep GPU-related rows when GPU was requested, even if runtime fell back to CPU.
            rows = [
                ["Requested mode", _RUNTIME["requested_mode"]],
                ["Requested GPU backend", _RUNTIME["requested_gpu_backend"]],
                ["Active mode", _RUNTIME["active_mode"]],
                ["Active GPU backend", _RUNTIME["active_gpu_backend"]],
                ["Provider", _RUNTIME["provider"]],
                ["Reason", _RUNTIME["reason"]],
                ["GPU batch lines", _RUNTIME["gpu_batch_lines"]],
                ["GPU batch grid", _RUNTIME["gpu_batch_grid"]],
            ]
        print()
        print(tabulate(rows, headers=["Field", "Value"], tablefmt="fancy_grid"))
        print()

    return dict(_RUNTIME)


def requested_mode() -> str:
    return str(_RUNTIME['requested_mode'])


def active_mode() -> str:
    return str(_RUNTIME['active_mode'])


def requested_gpu_backend() -> str:
    return str(_RUNTIME['requested_gpu_backend'])


def active_gpu_backend() -> str:
    return str(_RUNTIME['active_gpu_backend'])


def get_provider() -> Optional[str]:
    return _RUNTIME['provider']


def using_gpu() -> bool:
    return active_mode() == 'GPU' and get_provider() is not None


def is_cuda_active() -> bool:
    return using_gpu() and active_gpu_backend() == 'CUDA'


def is_mps_active() -> bool:
    return using_gpu() and active_gpu_backend() == 'MPS'


def get_cupy() -> Optional[Any]:
    if _RUNTIME['cupy'] is None:
        _try_import_cupy()
    return _RUNTIME['cupy']


def get_torch() -> Optional[Any]:
    if _RUNTIME['torch'] is None:
        _try_import_torch()
    return _RUNTIME['torch']


def get_torch_device() -> Optional[str]:
    return _RUNTIME['torch_device']


def get_gpu_batch_lines() -> int:
    return int(_RUNTIME['gpu_batch_lines'])


def get_gpu_batch_grid() -> int:
    return int(_RUNTIME['gpu_batch_grid'])


def to_numpy(arr: Any) -> np.ndarray:
    provider = get_provider()
    if provider == 'cupy':
        cp = get_cupy()
        if cp is not None and isinstance(arr, cp.ndarray):
            return cp.asnumpy(arr)
    elif provider == 'torch':
        torch = get_torch()
        if torch is not None and isinstance(arr, torch.Tensor):
            return arr.detach().cpu().numpy()
    return np.asarray(arr)


def free_gpu_memory() -> None:
    if not using_gpu():
        return

    provider = get_provider()
    if provider == 'cupy':
        cp = get_cupy()
        if cp is None:
            return
        try:
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass
    elif provider == 'torch':
        torch = get_torch()
        if torch is None:
            return
        try:
            if is_cuda_active() and torch.cuda.is_available():
                torch.cuda.empty_cache()
            elif is_mps_active():
                mps_mod = getattr(torch, 'mps', None)
                if mps_mod is not None and hasattr(mps_mod, 'empty_cache'):
                    mps_mod.empty_cache()
        except Exception:
            pass


def _torch_dtype(torch, device: str):
    if str(device).startswith('mps'):
        return torch.float32
    return torch.float64


def _cupy_arr(x, cp):
    return cp.asarray(np.asarray(x), dtype=cp.float64)


def _torch_arr(x, torch, device):
    return torch.as_tensor(np.asarray(x), dtype=_torch_dtype(torch, device), device=device)


def _backend_arrays(*vals) -> Optional[Tuple[str, Any, list]]:
    if not using_gpu():
        return None
    provider = get_provider()
    if provider == 'cupy':
        cp = get_cupy()
        if cp is None:
            return None
        return ('cupy', cp, [_cupy_arr(v, cp) for v in vals])
    if provider == 'torch':
        torch = get_torch()
        device = get_torch_device()
        if torch is None or device is None:
            return None
        return ('torch', torch, [_torch_arr(v, torch, device) for v in vals])
    return None


def _to_numpy(provider: str, mod, value):
    if provider == 'cupy':
        return mod.asnumpy(value)
    return value.detach().cpu().numpy()
