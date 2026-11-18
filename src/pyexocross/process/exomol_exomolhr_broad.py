"""Shared ExoMol and ExoMolHR pressure-broadening recipes."""
import glob
import os

import numpy as np
import pandas as pd


BROAD_DEFAULT_GAMMA_L = 0.07
BROAD_DEFAULT_N_AIR = 0.5
RECIPE_QNS = {
    'a0': ('J"',),
    'a1': ('J"', 'K"'),
    'm0': ('m',),
    'm1': ('m', 'K"'),
    'a5': ("J'", "Ka'", "Kc'", 'J"', 'Ka"', 'Kc"'),
}
RECIPE_LEVELS = ((5, ('a5',)), (1, ('a1', 'm1')), (0, ('a0', 'm0')))


def broad_recipe_qn_columns(code):
    """Return the quantum numbers stored after code, gamma and n."""
    return RECIPE_QNS.get(str(code).lower().strip()[:2], ())


def broad_required_line_columns(broad_dfs, available_columns=None):
    """Select line columns needed by the present broadening recipes."""
    available = set(available_columns) if available_columns is not None else None
    recipes = set()
    for table in broad_dfs:
        if 'code' in table:
            recipes.update(table['code'].astype(str).str.lower().str.strip().str[:2])
        else:
            recipes.add('a0')
    required = []
    for recipe in sorted(recipes):
        for name in RECIPE_QNS.get(recipe, ()):
            if name == 'm':
                required.extend(col for col in ("J'", 'J"') if available is None or col in available)
            elif available is None or name in available:
                required.append(name)
    return list(dict.fromkeys(required))


def read_broad_file(fname):
    """Read a .broad table, retaining every recipe-dependent QN column."""
    with open(fname) as handle:
        widths = [len(line.split()) for line in handle if line.strip() and not line.lstrip().startswith('#')]
    if not widths or min(widths) < 4:
        raise ValueError(f'Invalid broadening file {fname}: expected at least 4 columns.')
    names = ['code', 'gamma_L', 'n_air'] + [f'q{i}' for i in range(1, max(widths) - 2)]
    table = pd.read_csv(fname, sep=r'\s+', header=None, names=names, comment='#')
    table['code'] = table['code'].astype(str)
    for col in names[1:]:
        table[col] = pd.to_numeric(table[col], errors='coerce')
    return table


def read_broad(read_path):
    """Load one native .broad table or default per active broadener."""
    from pyexocross.core import data_info, broadeners, ratios, gamma_broadening, n_broadening

    broad, ratio, tables = [], [], []
    for i, amount in enumerate(ratios):
        if amount == 0:
            continue
        name = str(broadeners[i])
        if gamma_broadening or name.upper().startswith('DEF'):
            gamma = gamma_broadening[i] if gamma_broadening else BROAD_DEFAULT_GAMMA_L
            exponent = n_broadening[i] if gamma_broadening else BROAD_DEFAULT_N_AIR
            table = pd.DataFrame([{'code': 'Default', 'gamma_L': gamma, 'n_air': exponent}])
            broad.append(name if gamma_broadening else 'Default')
        else:
            pattern = os.path.join(read_path, data_info[0], '**', '*.broad')
            expected = f'{data_info[1]}__{name}.broad'.lower()
            files = sorted(
                path for path in glob.glob(pattern, recursive=True)
                if os.path.basename(path).lower() == expected
            )
            if len(files) != 1:
                raise ValueError(f'Expected one {name} broadening file, found {len(files)}.')
            table = read_broad_file(files[0])
            broad.append(name)
        ratio.append(amount)
        tables.append(table)
    return broad, ratio, len(broad), tables


def qn_number(values):
    return pd.to_numeric(values, errors='coerce').to_numpy(dtype=float)


def line_m(lines):
    """Calculate signed m from upper and lower J; other branches stay unmatched."""
    jpp = qn_number(lines['J"'])
    dj = qn_number(lines["J'"]) - jpp
    m = np.full(len(lines), np.nan)
    m[np.isclose(dj, -1)] = -jpp[np.isclose(dj, -1)]
    m[np.isclose(dj, 0)] = jpp[np.isclose(dj, 0)]
    m[np.isclose(dj, 1)] = (jpp + 1)[np.isclose(dj, 1)]
    return m


def recipe_values(recipe, table, lines):
    """Return gamma and n matched by one recipe, or None if QNs are absent."""
    names = RECIPE_QNS[recipe]
    line_columns = [col for col in names if col != 'm']
    if 'm' in names:
        line_columns.extend(("J'", 'J"'))
    file_columns = [f'q{i}' for i in range(1, len(names) + 1)]
    if not set(line_columns).issubset(lines.columns) or not set(file_columns).issubset(table.columns):
        return None

    source = table.dropna(subset=file_columns + ['gamma_L', 'n_air']).copy()
    if source.empty:
        return None
    keys = [qn_number(source[col]) for col in file_columns]
    target = [line_m(lines) if name == 'm' else qn_number(lines[name]) for name in names]
    if 'm' in names and not np.any(keys[0] < 0):
        target[0] = np.abs(target[0])
    if recipe == 'a0':
        target[0] = np.minimum(target[0], np.max(keys[0]))

    key_names = [f'key{i}' for i in range(len(names))]
    for key_name, values in zip(key_names, keys):
        source[key_name] = values
    # Repeated keys in existing .broad files use the final row.
    source = source.drop_duplicates(key_names, keep='last')
    lookup = source.set_index(key_names)[['gamma_L', 'n_air']]
    if len(names) == 1:
        index = pd.Index(target[0], name=key_names[0])
    else:
        index = pd.MultiIndex.from_arrays(target, names=key_names)
    result = lookup.reindex(index)
    return result['gamma_L'].to_numpy(), result['n_air'].to_numpy()


def extract_broad(broad_df, st_df):
    """Match ExoMol diet rows to transitions, with numeric recipe priority."""
    if len(st_df) == 0:
        empty = pd.Series([], index=st_df.index, dtype=float)
        return empty, empty.copy()
    table = broad_df.copy()
    if 'q1' not in table and 'Jpp' in table:
        table['q1'] = table['Jpp']
    if 'code' not in table:
        table['code'] = 'a0'
    table['recipe'] = table['code'].astype(str).str.lower().str.strip().str[:2]
    for col in ['gamma_L', 'n_air'] + [name for name in table if name.startswith('q')]:
        table[col] = pd.to_numeric(table[col], errors='coerce')

    gamma = np.full(len(st_df), BROAD_DEFAULT_GAMMA_L)
    exponent = np.full(len(st_df), BROAD_DEFAULT_N_AIR)
    unmatched = np.ones(len(st_df), dtype=bool)
    for _, recipes in RECIPE_LEVELS:
        tier_gamma = np.full(len(st_df), np.nan)
        tier_exponent = np.full(len(st_df), np.nan)
        tier_recipe = np.full(len(st_df), '', dtype=object)
        for recipe in recipes:
            rows = table[table['recipe'] == recipe]
            if rows.empty:
                continue
            result = recipe_values(recipe, rows, st_df)
            if result is None:
                continue
            current_gamma, current_exponent = result
            valid = np.isfinite(current_gamma) & np.isfinite(current_exponent) & unmatched
            overlap = valid & np.isfinite(tier_gamma)
            conflict = overlap & (
                ~np.isclose(tier_gamma, current_gamma, equal_nan=True)
                | ~np.isclose(tier_exponent, current_exponent, equal_nan=True)
            )
            if conflict.any():
                index = st_df.index[np.flatnonzero(conflict)[0]]
                raise ValueError(
                    f'Conflicting {tier_recipe[np.flatnonzero(conflict)[0]]} and {recipe} '
                    f'broadening rows for transition {index}.'
                )
            fill = valid & ~np.isfinite(tier_gamma)
            tier_gamma[fill] = current_gamma[fill]
            tier_exponent[fill] = current_exponent[fill]
            tier_recipe[fill] = recipe
        fill = unmatched & np.isfinite(tier_gamma)
        gamma[fill] = tier_gamma[fill]
        exponent[fill] = tier_exponent[fill]
        unmatched[fill] = False

    return pd.Series(gamma, index=st_df.index), pd.Series(exponent, index=st_df.index)
