"""GPU cross-section kernels for CUDA / MPS backends."""

from __future__ import annotations

from typing import Any, Optional

import numpy as np
from scipy.special import roots_hermite

from pyexocross.gpu.base_gpu import (
    free_gpu_memory,
    get_cupy,
    get_gpu_batch_grid,
    get_gpu_batch_lines,
    get_provider,
    get_torch,
    get_torch_device,
    using_gpu,
)

_GPU_NOTICE_KEYS = set()
_EPS = float(np.finfo(np.float64).eps)


def _notify_gpu_fallback(kind, message, key_type='fallback'):
    key = (key_type, kind)
    if key in _GPU_NOTICE_KEYS:
        return
    _GPU_NOTICE_KEYS.add(key)
    print(f"[GPU fallback] {kind}: {message}")


def _window_bounds(wn_grid, v, cutoff):
    if cutoff == 'None':
        start = max(0, wn_grid.searchsorted(np.min(v)) - 1)
        end = min(wn_grid.searchsorted(np.max(v)), len(wn_grid))
    else:
        wing = float(cutoff)
        start = max(0, wn_grid.searchsorted(np.min(v) - wing) - 1)
        end = min(wn_grid.searchsorted(np.max(v) + wing), len(wn_grid))
    return start, end


def _torch_float_dtype(torch, device):
    if str(device).startswith('mps'):
        return torch.float32
    return torch.float64


def _torch_complex_dtype(torch, device):
    if str(device).startswith('mps'):
        return torch.complex64
    return torch.complex128


def _finite_or_zero(provider, mod, arr):
    return mod.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def _safe_positive(provider, mod, arr, eps=_EPS):
    return mod.where(arr > eps, arr, mod.full_like(arr, eps))


def _safe_nonzero(provider, mod, arr, eps=_EPS):
    sign = mod.where(arr < 0, -1.0, 1.0)
    return mod.where(mod.abs(arr) > eps, arr, sign * eps)


def _cupy_humlicek_profile(cp, dv, alpha, gamma):
    alpha = cp.where(alpha > _EPS, alpha, cp.full_like(alpha, _EPS))
    gamma = cp.where(gamma > _EPS, gamma, cp.full_like(gamma, _EPS))
    sqrtln2 = np.sqrt(np.log(2.0))
    inv_sqrt_pi = 1.0 / np.sqrt(np.pi)
    sqrtln2_inv_pi = np.sqrt(np.log(2.0) / np.pi)

    x = dv * sqrtln2 / alpha[None, :]
    y = gamma[None, :] * sqrtln2 / alpha[None, :]
    t = y - 1j * x
    s = cp.abs(x) + y
    u = t**2
    ybound = 0.195 * cp.abs(x) - 0.176

    w = cp.zeros_like(t, dtype=cp.complex128)

    mask1 = s >= 15
    w1 = t * inv_sqrt_pi / (0.5 + t**2)
    w = cp.where(mask1, w1, w)

    mask2 = (s >= 5.5) & (s < 15)
    w2 = (t * (1.4104739589 + u * inv_sqrt_pi)) / (0.75 + u * (3.0 + u))
    w = cp.where(mask2, w2, w)

    mask3 = (s < 5.5) & (y >= ybound)
    w3 = (
        16.4955
        + t * (20.20933 + t * (11.96482 + t * (3.778987 + 0.5642236 * t)))
    ) / (
        16.4955
        + t * (38.82363 + t * (39.27121 + t * (21.69274 + t * (6.699398 + t))))
    )
    w = cp.where(mask3, w3, w)

    mask4 = (s < 5.5) & (y < ybound)
    nom = t * (
        36183.31
        - u * (3321.99 - u * (1540.787 - u * (219.031 - u * (35.7668 - u * (1.320522 - u * 0.56419)))))
    )
    den = (
        32066.6
        - u * (24322.8 - u * (9022.23 - u * (2186.18 - u * (364.219 - u * (61.5704 - u * (1.84144 - u))))))
    )
    w4 = cp.exp(u) - nom / den
    w = cp.where(mask4, w4, w)

    out = cp.real(w) / alpha[None, :] * sqrtln2_inv_pi
    return cp.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


def _torch_humlicek_profile(torch, dv, alpha, gamma):
    device = dv.device.type
    fdtype = _torch_float_dtype(torch, device)
    alpha = torch.where(alpha > _EPS, alpha, torch.full_like(alpha, _EPS))
    gamma = torch.where(gamma > _EPS, gamma, torch.full_like(gamma, _EPS))
    sqrtln2 = float(np.sqrt(np.log(2.0)))
    sqrtln2_inv_pi = float(np.sqrt(np.log(2.0) / np.pi))

    try:
        cdtype = _torch_complex_dtype(torch, device)
        inv_sqrt_pi = float(1.0 / np.sqrt(np.pi))

        x = dv * sqrtln2 / alpha[None, :]
        y = gamma[None, :] * sqrtln2 / alpha[None, :]
        t = y.to(cdtype) - (1j) * x.to(cdtype)
        s = torch.abs(x) + y
        u = t**2
        ybound = 0.195 * torch.abs(x) - 0.176

        w = torch.zeros_like(t, dtype=cdtype)

        mask1 = s >= 15
        w1 = t * inv_sqrt_pi / (0.5 + t**2)
        w = torch.where(mask1, w1, w)

        mask2 = (s >= 5.5) & (s < 15)
        w2 = (t * (1.4104739589 + u * inv_sqrt_pi)) / (0.75 + u * (3.0 + u))
        w = torch.where(mask2, w2, w)

        mask3 = (s < 5.5) & (y >= ybound)
        w3 = (
            16.4955
            + t * (20.20933 + t * (11.96482 + t * (3.778987 + 0.5642236 * t)))
        ) / (
            16.4955
            + t * (38.82363 + t * (39.27121 + t * (21.69274 + t * (6.699398 + t))))
        )
        w = torch.where(mask3, w3, w)

        mask4 = (s < 5.5) & (y < ybound)
        nom = t * (
            36183.31
            - u * (3321.99 - u * (1540.787 - u * (219.031 - u * (35.7668 - u * (1.320522 - u * 0.56419)))))
        )
        den = (
            32066.6
            - u * (24322.8 - u * (9022.23 - u * (2186.18 - u * (364.219 - u * (61.5704 - u * (1.84144 - u))))))
        )
        w4 = torch.exp(u) - nom / den
        w = torch.where(mask4, w4, w)

        out = torch.real(w) / alpha[None, :] * sqrtln2_inv_pi
        out = torch.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
        return out.to(fdtype)
    except Exception:
        # Complex math may be unavailable on some MPS builds; use pseudo-Voigt
        # on-GPU approximation to keep GPU execution active.
        gaussian = sqrtln2_inv_pi / alpha[None, :] * torch.exp(-np.log(2.0) * (dv / alpha[None, :]) ** 2)
        lorentzian = gamma[None, :] / np.pi / (dv**2 + gamma[None, :] ** 2)
        hV = (
            alpha**5
            + 2.69269 * alpha**4 * gamma
            + 2.42843 * alpha**3 * gamma**2
            + 4.47163 * alpha**2 * gamma**3
            + 0.07842 * alpha * gamma**4
            + gamma**5
        ) ** 0.2
        eta = 1.36603 * (gamma / hV) - 0.47719 * (gamma / hV) ** 2 + 0.11116 * (gamma / hV) ** 3
        out = eta[None, :] * lorentzian + (1.0 - eta[None, :]) * gaussian
        out = torch.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
        return out.to(fdtype)


def _cupy_voigt_profile(cp, dv, sigma, gamma):
    sigma = cp.where(sigma > _EPS, sigma, cp.full_like(sigma, _EPS))
    gamma = cp.where(gamma > _EPS, gamma, cp.full_like(gamma, _EPS))
    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    inv_sqrt2pi = 1.0 / np.sqrt(2.0 * np.pi)
    try:
        from cupyx.scipy.special import wofz as cp_wofz  # type: ignore
        z = (dv + 1j * gamma[None, :]) / sigma[None, :] * inv_sqrt2
        wz = cp_wofz(z)
        out = cp.real(wz) / sigma[None, :] * inv_sqrt2pi
        return cp.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
    except Exception:
        alpha = sigma * np.sqrt(2.0 * np.log(2.0))
        return _cupy_humlicek_profile(cp, dv, alpha, gamma)


def _torch_voigt_profile(torch, dv, sigma, gamma):
    # Torch currently has no native wofz on all devices; use Humlicek approximation.
    sigma = torch.where(sigma > _EPS, sigma, torch.full_like(sigma, _EPS))
    gamma = torch.where(gamma > _EPS, gamma, torch.full_like(gamma, _EPS))
    alpha = sigma * np.sqrt(2.0 * np.log(2.0))
    return _torch_humlicek_profile(torch, dv, alpha, gamma)


def _to_backend_array(provider, mod, x, *, is_grid=False):
    arr = np.asarray(x)
    if provider == 'cupy':
        return mod.asarray(arr, dtype=mod.float64)

    # torch
    device = get_torch_device()
    dtype = _torch_float_dtype(mod, device)
    return mod.as_tensor(arr, dtype=dtype, device=device)


def gpu_cross_section(
    kind,
    wn_grid,
    v,
    coef,
    cutoff,
    *,
    alpha=None,
    gamma=None,
    sigma=None,
    eta=None,
    bin_size2=2.0,
    inv_bin_size_pi_half=1.0,
):
    """Compute cross section on GPU. Returns None when GPU is inactive."""
    if not using_gpu():
        return None

    supported = {
        'doppler',
        'lorentzian',
        'scipy_voigt',
        'scipy_wofz_voigt',
        'humlicek_voigt',
        'pseudo_voigt',
        'binned_gaussian',
        'binned_lorentzian',
        'binned_voigt',
    }
    if kind not in supported:
        _notify_gpu_fallback(kind, 'profile currently uses CPU implementation.', key_type='unsupported')
        return None

    provider = get_provider()
    if provider == 'cupy':
        mod = get_cupy()
    elif provider == 'torch':
        mod = get_torch()
    else:
        return None

    if mod is None:
        return None

    try:
        wn_grid_np = np.asarray(wn_grid, dtype=np.float64)
        v_np = np.asarray(v, dtype=np.float64)
        coef_np = np.asarray(coef, dtype=np.float64)
        if len(v_np) == 0:
            return np.zeros_like(wn_grid_np)

        alpha_np = np.asarray(alpha, dtype=np.float64) if alpha is not None else None
        gamma_np = np.asarray(gamma, dtype=np.float64) if gamma is not None else None
        sigma_np = np.asarray(sigma, dtype=np.float64) if sigma is not None else None
        eta_np = np.asarray(eta, dtype=np.float64) if eta is not None else None

        # Keep line filtering consistent with CPU behavior: only drop clearly
        # invalid line centers/coefficients. Width parameters are sanitized
        # below instead of filtering out whole lines.
        valid = np.isfinite(v_np) & np.isfinite(coef_np)
        if not np.any(valid):
            return np.zeros_like(wn_grid_np)

        v_np = v_np[valid]
        coef_np = np.nan_to_num(coef_np[valid], nan=0.0, posinf=0.0, neginf=0.0)
        if alpha_np is not None:
            alpha_np = np.nan_to_num(alpha_np[valid], nan=_EPS, posinf=_EPS, neginf=_EPS)
            alpha_np = np.where(alpha_np > _EPS, alpha_np, _EPS)
        if gamma_np is not None:
            gamma_np = np.nan_to_num(gamma_np[valid], nan=_EPS, posinf=_EPS, neginf=_EPS)
            gamma_np = np.where(gamma_np > _EPS, gamma_np, _EPS)
        if sigma_np is not None:
            sigma_np = np.nan_to_num(sigma_np[valid], nan=_EPS, posinf=_EPS, neginf=_EPS)
            sigma_np = np.where(sigma_np > _EPS, sigma_np, _EPS)
        if eta_np is not None:
            eta_np = np.nan_to_num(eta_np[valid], nan=0.5, posinf=0.5, neginf=0.5)

        start, end = _window_bounds(wn_grid_np, v_np, cutoff)
        xsec = np.zeros_like(wn_grid_np, dtype=np.float64)
        if start >= end:
            return xsec

        line_batch = max(1, min(get_gpu_batch_lines(), len(v_np)))
        grid_batch = max(1, get_gpu_batch_grid())
        cutoff_val = None if cutoff == 'None' else float(cutoff)

        sqrtln2 = np.sqrt(np.log(2.0))
        sqrtln2_inv_pi = np.sqrt(np.log(2.0) / np.pi)
        negln2 = -np.log(2.0)
        pi_val = np.pi
        bin_half = float(bin_size2) / 2.0

        # global normalization values for binned line profiles
        bin_size = float(bin_size2) / 2.0
        wngrid_start = wn_grid_np[start]
        wngrid_end = wn_grid_np[end - 1]

        if kind == 'binned_lorentzian':
            if gamma_np is None:
                return None
            denom = np.arctan((wngrid_end - v_np) / gamma_np) - np.arctan((wngrid_start - v_np) / gamma_np)
            safe_denom = np.where(np.abs(denom) > np.finfo(float).eps, denom, np.finfo(float).eps)
            bnorm_np = 1.0 / safe_denom / bin_size
        else:
            bnorm_np = None

        if kind == 'binned_voigt':
            if sigma_np is None or gamma_np is None:
                return None
            nquad = 20
            roots, weights = roots_hermite(nquad, mu=False)
        else:
            roots, weights = None, None

        for g0 in range(start, end, grid_batch):
            g1 = min(g0 + grid_batch, end)
            wn_block = _to_backend_array(provider, mod, wn_grid_np[g0:g1])[:, None]
            if provider == 'torch':
                block_sum = mod.zeros((g1 - g0,), dtype=wn_block.dtype, device=wn_block.device)
            else:
                block_sum = mod.zeros((g1 - g0,), dtype=wn_block.dtype)

            for l0 in range(0, len(v_np), line_batch):
                l1 = min(l0 + line_batch, len(v_np))
                v_g = _to_backend_array(provider, mod, v_np[l0:l1])
                coef_g = _to_backend_array(provider, mod, coef_np[l0:l1])
                dv = wn_block - v_g[None, :]

                alpha_g = _to_backend_array(provider, mod, alpha_np[l0:l1]) if alpha_np is not None else None
                gamma_g = _to_backend_array(provider, mod, gamma_np[l0:l1]) if gamma_np is not None else None
                sigma_g = _to_backend_array(provider, mod, sigma_np[l0:l1]) if sigma_np is not None else None
                eta_g = _to_backend_array(provider, mod, eta_np[l0:l1]) if eta_np is not None else None

                if kind == 'doppler':
                    alpha_safe = _safe_positive(provider, mod, alpha_g)
                    profile = sqrtln2_inv_pi / alpha_safe[None, :] * mod.exp(negln2 * (dv / alpha_safe[None, :]) ** 2)
                elif kind == 'lorentzian':
                    gamma_safe = _safe_positive(provider, mod, gamma_g)
                    profile = gamma_safe[None, :] / pi_val / (dv**2 + gamma_safe[None, :] ** 2)
                elif kind in ('scipy_voigt', 'scipy_wofz_voigt'):
                    if provider == 'cupy':
                        profile = _cupy_voigt_profile(mod, dv, sigma_g, gamma_g)
                    else:
                        profile = _torch_voigt_profile(mod, dv, sigma_g, gamma_g)
                elif kind == 'humlicek_voigt':
                    if provider == 'cupy':
                        profile = _cupy_humlicek_profile(mod, dv, alpha_g, gamma_g)
                    else:
                        profile = _torch_humlicek_profile(mod, dv, alpha_g, gamma_g)
                elif kind == 'pseudo_voigt':
                    alpha_safe = _safe_positive(provider, mod, alpha_g)
                    gamma_safe = _safe_positive(provider, mod, gamma_g)
                    gaussian = sqrtln2_inv_pi / alpha_safe[None, :] * mod.exp(negln2 * (dv / alpha_safe[None, :]) ** 2)
                    lorentzian = gamma_safe[None, :] / pi_val / (dv**2 + gamma_safe[None, :] ** 2)
                    profile = eta_g[None, :] * lorentzian + (1.0 - eta_g[None, :]) * gaussian
                elif kind == 'binned_gaussian':
                    alpha_safe = _safe_positive(provider, mod, alpha_g)
                    if provider == 'cupy':
                        try:
                            erfxpos = mod.erf(sqrtln2 * (dv + bin_half) / alpha_safe[None, :])
                            erfxneg = mod.erf(sqrtln2 * (dv - bin_half) / alpha_safe[None, :])
                        except Exception:
                            from cupyx.scipy.special import erf as cp_erf  # type: ignore
                            erfxpos = cp_erf(sqrtln2 * (dv + bin_half) / alpha_safe[None, :])
                            erfxneg = cp_erf(sqrtln2 * (dv - bin_half) / alpha_safe[None, :])
                    else:
                        erfxpos = mod.special.erf(sqrtln2 * (dv + bin_half) / alpha_safe[None, :])
                        erfxneg = mod.special.erf(sqrtln2 * (dv - bin_half) / alpha_safe[None, :])
                    profile = erfxpos - erfxneg
                elif kind == 'binned_lorentzian':
                    gamma_safe = _safe_positive(provider, mod, gamma_g)
                    bnorm_g = _to_backend_array(provider, mod, bnorm_np[l0:l1])
                    profile = (
                        mod.atan((dv + bin_half) / gamma_safe[None, :])
                        - mod.atan((dv - bin_half) / gamma_safe[None, :])
                    ) * bnorm_g[None, :]
                elif kind == 'binned_voigt':
                    sigma_safe = _safe_positive(provider, mod, sigma_g)
                    gamma_safe = _safe_positive(provider, mod, gamma_g)
                    # Build quadrature tensors: [g, l, q]
                    roots_g = _to_backend_array(provider, mod, roots)
                    weights_g = _to_backend_array(provider, mod, weights)

                    # bnormq: [l, q]
                    vxsigma = v_g[:, None] + roots_g[None, :] * sigma_safe[:, None]
                    bnormq_den = (
                        mod.atan((wngrid_end - vxsigma) / gamma_safe[:, None])
                        - mod.atan((wngrid_start - vxsigma) / gamma_safe[:, None])
                    )
                    bnormq = 1.0 / _safe_nonzero(provider, mod, bnormq_den)

                    dvx = dv[:, :, None] - roots_g[None, None, :] * sigma_safe[None, :, None]
                    lorenz = (
                        mod.atan((dvx + bin_half) / gamma_safe[None, :, None])
                        - mod.atan((dvx - bin_half) / gamma_safe[None, :, None])
                    )
                    profile = mod.sum(weights_g[None, None, :] * bnormq[None, :, :] * lorenz, axis=2)
                else:
                    return None

                profile = _finite_or_zero(provider, mod, profile)
                if cutoff_val is not None:
                    mask = mod.abs(dv) <= cutoff_val
                    profile = mod.where(mask, profile, 0.0)
                    profile = _finite_or_zero(provider, mod, profile)

                contribution = _finite_or_zero(provider, mod, profile * coef_g[None, :])
                block_sum += mod.sum(contribution, axis=1)
                block_sum = _finite_or_zero(provider, mod, block_sum)

            if provider == 'cupy':
                xsec[g0:g1] = mod.asnumpy(_finite_or_zero(provider, mod, block_sum))
            else:
                xsec[g0:g1] = _finite_or_zero(provider, mod, block_sum).detach().cpu().numpy()

        if kind == 'binned_gaussian':
            xsec[start:end] /= float(bin_size2)
        elif kind == 'binned_voigt':
            xsec[start:end] *= float(inv_bin_size_pi_half)
        xsec = np.nan_to_num(xsec, nan=0.0, posinf=0.0, neginf=0.0)

        free_gpu_memory()
        return xsec

    except Exception as exc:
        _notify_gpu_fallback(kind, f"GPU path failed ({exc}); using CPU implementation.", key_type='error')
        free_gpu_memory()
        return None
