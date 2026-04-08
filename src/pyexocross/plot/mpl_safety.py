"""Matplotlib safety helpers for very large line/path plots."""

from __future__ import annotations

import matplotlib as mpl


def configure_mpl_large_path_rendering() -> None:
    """Tune Agg/path settings to avoid OverflowError on large curves."""
    mpl.rcParams["agg.path.chunksize"] = 20000
    mpl.rcParams["path.simplify"] = True
    mpl.rcParams["path.simplify_threshold"] = 1.0
