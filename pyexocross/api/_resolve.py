"""
Database metadata resolver for PyExoCross API.

Auto-detects states file columns, format strings, and database-specific
metadata (species_main_id, species_sub_id, abundance, mass, capability flags)
by reading actual database definition files.

This allows the Python API to work without an .inp file by deriving the
same metadata that inp_para() reads from .inp + database files.
"""
import os
import re
import numpy as np
import pandas as pd


def resolve_database_metadata(database, read_path, data_info, species_id=None):
    """
    Auto-detect database-specific metadata from definition files.

    Reads the ExoMol/ExoAtom definition file or HITRAN molparam table to
    determine column names/formats, capability flags, and physical constants.

    Parameters
    ----------
    database : str
        Database name: 'ExoMol', 'ExoAtom', 'HITRAN', or 'HITEMP'.
    read_path : str
        Path to the database root directory (ExoMol/ExoAtom) or the .par
        file (HITRAN/HITEMP).
    data_info : list of str
        [molecule, isotopologue, dataset] for ExoMol/HITRAN/HITEMP, or
        [atom, dataset] for ExoAtom.
    species_id : int, optional
        HITRAN molecule-isotopologue ID (e.g. 81 for NO 14N-16O).
        Required for HITRAN/HITEMP databases.

    Returns
    -------
    dict
        Dictionary with keys:
        - states_col : list of str
        - states_fmt : list of str
        - check_uncertainty : bool
        - check_lifetime : bool
        - check_gfactor : bool
        - check_predissoc : bool
        - species_main_id : int
        - species_sub_id : int
        - abundance : float
        - mass : float

    Raises
    ------
    ValueError
        If the database name is not recognized.
    FileNotFoundError
        If required definition or metadata files cannot be found.
    """
    if database == 'ExoMol':
        return _resolve_exomol(read_path, data_info, species_id)
    elif database == 'ExoAtom':
        return _resolve_exoatom(read_path, data_info, species_id)
    elif database in ('HITRAN', 'HITEMP'):
        return _resolve_hitran(read_path, data_info, species_id)
    else:
        raise ValueError(
            f"Unsupported database '{database}'. "
            "Choose from 'ExoMol', 'ExoAtom', 'HITRAN', or 'HITEMP'."
        )


def _resolve_exomol(read_path, data_info, species_id=None):
    """Resolve metadata from ExoMol definition files (.def.json or .def)."""
    base_dir = os.path.join(read_path, *data_info)
    base_name = '__'.join(data_info[-2:])

    # Give default values, then override with species_id
    species_main_id     = 0
    species_sub_id = 0
    if species_id is not None:
        species_main_id = int(species_id / 10)
        species_sub_id = species_id - species_main_id * 10

    # Try JSON definition file first
    json_path = os.path.join(base_dir, base_name + '.def.json')
    if os.path.exists(json_path):
        meta = _parse_exomol_json_def(json_path)
        meta.setdefault('species_main_id', species_main_id)
        meta.setdefault('species_sub_id', species_sub_id)
        return meta

    # Fall back to plain-text definition file
    def_path = os.path.join(base_dir, base_name + '.def')
    if os.path.exists(def_path):
        meta = _parse_exomol_text_def(def_path)
        meta.setdefault('species_main_id', species_main_id)
        meta.setdefault('species_sub_id', species_sub_id)
        return meta

    raise FileNotFoundError(
        f"ExoMol definition file not found. "
        f"Looked for: {json_path} and {def_path}"
    )


def _parse_exomol_json_def(json_path):
    """Parse ExoMol JSON definition file."""
    def_df = pd.read_json(json_path, orient='columns')
    states_fields = def_df['dataset']['states']['states_file_fields']
    states_col = [f['name'] for f in states_fields]
    states_fmt = [f['cfmt'] for f in states_fields]

    check_uncertainty = _safe_get(def_df, ['dataset', 'states', 'uncertainties_available'], False)
    check_predissoc = _safe_get(def_df, ['dataset', 'predis'], False) or False
    check_lifetime = _safe_get(def_df, ['dataset', 'states', 'lifetime_available'], False)
    check_gfactor = _safe_get(def_df, ['dataset', 'states', 'lande_g_available'], False)

    abundance = 1.0
    mass = def_df['isotopologue']['mass_in_Da']

    return {
        'states_col': states_col,
        'states_fmt': states_fmt,
        'check_uncertainty': bool(check_uncertainty),
        'check_lifetime': bool(check_lifetime),
        'check_gfactor': bool(check_gfactor),
        'check_predissoc': bool(check_predissoc),
        'abundance': abundance,
        'mass': float(mass),
    }


def _parse_exomol_text_def(def_path):
    """Parse ExoMol plain-text definition file."""
    def_df = pd.read_csv(
        def_path, sep='\\s+', usecols=[0, 1, 2, 3, 4],
        names=['0', '1', '2', '3', '4'], header=None,
    )
    states_cfmt_df = def_df[def_df['3'].isin(['Format'])]
    states_col = def_df.iloc[states_cfmt_df.index - 1]['0'].tolist()
    states_fmt = states_cfmt_df['1'].tolist()

    def _bool_field(col_name, default=False):
        try:
            return bool(int(def_df[def_df['2'].isin([col_name])]['0'].values[0]))
        except (IndexError, ValueError):
            return default

    check_uncertainty = _bool_field('Uncertainty')
    check_predissoc = _bool_field('Predissociative')
    check_lifetime = _bool_field('Lifetime')

    # g-factor uses column '3' instead of '2'
    try:
        check_gfactor = bool(int(def_df[def_df['3'].isin(['g-factor'])]['0'].values[0]))
    except (IndexError, ValueError):
        check_gfactor = False

    abundance = 1.0
    mass = float(def_df[def_df['4'].isin(['mass'])]['0'].values[0])

    return {
        'states_col': states_col,
        'states_fmt': states_fmt,
        'check_uncertainty': check_uncertainty,
        'check_lifetime': check_lifetime,
        'check_gfactor': check_gfactor,
        'check_predissoc': check_predissoc,
        'abundance': abundance,
        'mass': mass,
    }


def _resolve_exoatom(read_path, data_info, species_id=None):
    """Resolve metadata from ExoAtom definition files (.adef.json)."""
    base_dir = os.path.join(read_path, *data_info)
    base_name = '__'.join(data_info[-2:])
    json_path = os.path.join(base_dir, base_name + '.adef.json')

    def_df = pd.read_json(json_path, orient='columns')
    states_fields = def_df['dataset']['states']['states_file_fields']
    states_col = [f['name'] for f in states_fields]
    states_fmt = [f['cfmt'] for f in states_fields]

    check_uncertainty = _safe_get(def_df, ['dataset', 'states', 'uncertainty_available'], False)
    check_predissoc = _safe_get(def_df, ['dataset', 'predis'], False) or False
    check_lifetime = _safe_get(def_df, ['dataset', 'states', 'lifetime_available'], False)
    check_gfactor = _safe_get(def_df, ['dataset', 'states', 'lande_g_available'], False)

    # Give default values, then override with species_id
    species_main_id     = 0
    species_sub_id = 0
    if species_id is not None:
        species_main_id = int(species_id / 10)
        species_sub_id = species_id - species_main_id * 10

    abundance = 1.0
    mass = float(def_df['species']['mass_in_Da'])

    return {
        'states_col': states_col,
        'states_fmt': states_fmt,
        'check_uncertainty': bool(check_uncertainty),
        'check_lifetime': bool(check_lifetime),
        'check_gfactor': bool(check_gfactor),
        'check_predissoc': bool(check_predissoc),
        'species_main_id': species_main_id,
        'species_sub_id': species_sub_id,
        'abundance': abundance,
        'mass': mass,
    }


def _resolve_hitran(read_path, data_info, species_id):
    """Resolve metadata from HITRAN molparam table."""
    if species_id is None:
        raise ValueError(
            "species_id is required for HITRAN/HITEMP databases. "
            "Provide the HITRAN molecule-isotopologue ID (e.g. 81 for NO 14N-16O)."
        )

    species_main_id = int(species_id / 10)
    species_sub_id = species_id - species_main_id * 10

    clean_path = read_path.rstrip('/')
    hitran_root = os.path.dirname(os.path.dirname(clean_path))
    molparam_filepath = os.path.join(hitran_root, 'molparam.txt')

    iso_meta_table = _parse_molparam(molparam_filepath, species_main_id)

    iso_meta_row = iso_meta_table[iso_meta_table['local ID'].isin([species_sub_id])]
    if iso_meta_row.empty:
        raise RuntimeError(
            f"Isotopologue ID {species_sub_id} not found for molecule ID {species_main_id} "
            f"in {molparam_filepath}"
        )

    abundance_val = iso_meta_row['Abundance'].iloc[0]
    if isinstance(abundance_val, str):
        abundance = float(abundance_val.replace('\xa0\u00d7\xa010', 'E'))
    else:
        abundance = float(abundance_val)
    mass = float(iso_meta_row['Molar Mass /g\u00b7mol-1'].iloc[0])

    return {
        'states_col': [],
        'states_fmt': [],
        'check_uncertainty': True,
        'check_lifetime': False,
        'check_gfactor': False,
        'check_predissoc': False,
        'species_main_id': species_main_id,
        'species_sub_id': species_sub_id,
        'abundance': abundance,
        'mass': mass,
    }


def _parse_molparam(molparam_filepath, species_main_id):
    """Parse HITRAN molparam.txt to extract isotopologue metadata."""
    local_ids = []
    abundances = []
    molar_masses = []
    current_mol_id = None

    try:
        with open(molparam_filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('Molecule #'):
                    continue
                m = re.search(r'\((\d+)\)\s*$', line)
                if m:
                    current_mol_id = int(m.group(1))
                    continue
                if current_mol_id != species_main_id:
                    continue
                parts = line.split()
                if len(parts) < 5:
                    continue
                try:
                    abundance_val = float(parts[1])
                    molar_mass_val = float(parts[4])
                except ValueError:
                    continue
                local_ids.append(len(local_ids) + 1)
                abundances.append(abundance_val)
                molar_masses.append(molar_mass_val)
    except FileNotFoundError:
        raise FileNotFoundError(
            f"HITRAN molparam.txt not found at: {molparam_filepath}. "
            f"Please ensure the file exists at the expected location."
        )

    if not local_ids:
        raise RuntimeError(
            f"No isotopologue entries found in molparam.txt for molecule ID {species_main_id}."
        )

    return pd.DataFrame({
        'local ID': np.array(local_ids, dtype=int),
        'Abundance': np.array(abundances, dtype=float),
        'Molar Mass /g\u00b7mol-1': np.array(molar_masses, dtype=float),
    })


def _safe_get(df, keys, default=False):
    """Safely traverse nested dictionary-like DataFrame access."""
    try:
        val = df
        for key in keys:
            val = val[key]
        return val
    except (KeyError, TypeError, IndexError):
        return default
