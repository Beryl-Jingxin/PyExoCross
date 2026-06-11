"""Download HITRAN line lists and partition functions."""

import os
import re
import requests
from bs4 import BeautifulSoup

from .common import default_url_path, download_url_to_path, resolve_save_path, write_urls


HITRAN_HOST = 'https://hitran.org'
HITRAN_ISO_META_URL = f'{HITRAN_HOST}/docs/iso-meta/'
HITRAN_PF_BASE_URL = f'{HITRAN_HOST}/data/Q'
_ISO_META_CACHE = None


def _canonical_formula(value):
    value = re.sub(r'[\(\)\-_\s]', '', str(value))
    value = re.sub(r'1H', 'H', value)
    return value


def _formula_cell_to_slug(formula_cell):
    """Convert a formula cell like ``14N16O`` to ``14N-16O``."""
    text = formula_cell.get_text('', strip=True)
    parts = re.findall(r'(\d*)([A-Z][a-z]?)(\d*)', text)
    slug_parts = []
    for isotope, element, count in parts:
        isotope = isotope or ('1' if element == 'H' else '')
        slug_parts.append(f'{isotope}{element}{count}')
    return '-'.join(slug_parts)


def _load_iso_meta(timeout=60):
    global _ISO_META_CACHE
    if _ISO_META_CACHE is not None:
        return _ISO_META_CACHE

    response = requests.get(HITRAN_ISO_META_URL, timeout=timeout)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, 'html.parser')

    iso_meta = {}
    for heading in soup.find_all('h4'):
        heading_text = heading.get_text(' ', strip=True)
        match = re.match(r'(\d+)\s*:\s*(.+)', heading_text)
        if match is None:
            continue
        molecule_id = int(match.group(1))
        molecule = match.group(2).replace(' ', '')
        table = heading.find_next('table')
        if table is None:
            continue
        for row in table.find_all('tr')[1:]:
            cells = row.find_all('td')
            if len(cells) < 8:
                continue
            global_id = int(cells[0].get_text(strip=True))
            local_id = int(cells[1].get_text(strip=True))
            isotopologue = _formula_cell_to_slug(cells[2])
            q_link = cells[7].find('a')
            q_href = q_link.get('href') if q_link is not None else f'/data/Q/q{global_id}.txt'
            q_url = q_href if q_href.startswith('http') else HITRAN_HOST + q_href
            iso_meta[(molecule, _canonical_formula(isotopologue))] = {
                'molecule': molecule,
                'molecule_id': molecule_id,
                'isotopologue': isotopologue,
                'local_id': local_id,
                'global_id': global_id,
                'q_url': q_url,
            }

    _ISO_META_CACHE = iso_meta
    return iso_meta


def _resolve_iso_meta(molecule, isotopologue, timeout=60):
    iso_meta = _load_iso_meta(timeout=timeout)
    key = (molecule, _canonical_formula(isotopologue))
    if key in iso_meta:
        return iso_meta[key]
    available = [
        row['isotopologue']
        for (mol, _), row in iso_meta.items()
        if mol == molecule
    ]
    raise ValueError(
        f'Cannot resolve HITRAN isotope metadata for {molecule}/{isotopologue}. '
        f'Available isotopologues: {available}'
    )


def _line_url(global_id, wn_range):
    wn_min, wn_max = wn_range
    return f'{HITRAN_HOST}/lbl/api?iso_ids_list={global_id}&numin={wn_min}&numax={wn_max}'


def _pf_url(global_id):
    return f'{HITRAN_PF_BASE_URL}/q{global_id}.txt'


def _target_base(save_path, molecule, isotopologue):
    return os.path.join(save_path, molecule, isotopologue, f'{molecule}__{isotopologue}')


def _normalize_hitran_configs(species_info, timeout=60):
    normalized = []
    for molecule, isotopologues in species_info.items():
        for isotopologue, config in isotopologues.items():
            if config is None:
                config = {}
            if not isinstance(config, dict):
                raise ValueError(f'HITRAN config for {molecule}/{isotopologue} must be a dict.')
            wn_range = config.get('wn_range')
            if wn_range is None:
                raise ValueError(f'HITRAN config for {molecule}/{isotopologue} requires wn_range.')
            meta = _resolve_iso_meta(molecule, isotopologue, timeout=timeout)
            normalized.append((molecule, isotopologue, meta, wn_range))
    return normalized


def get_hitran_urls(species_info, timeout=60):
    """Return HITRAN line API URLs and partition-function URLs."""
    urls = []
    configs = _normalize_hitran_configs(species_info, timeout=timeout)
    for _, _, meta, wn_range in configs:
        urls.extend([_line_url(meta['global_id'], wn_range), meta['q_url']])
    _print_hitran_urls(configs)
    return urls


def _print_hitran_urls(configs):
    """Print HITRAN URL groups in the same style as ExoMol downloads."""
    for molecule, isotopologue, meta, wn_range in configs:
        group_urls = [_line_url(meta['global_id'], wn_range), meta['q_url']]
        print(f'HITRAN - {molecule} - {isotopologue}: {len(group_urls)} file(s)')
        for url in group_urls:
            print(url)


def get_hitran_targets(species_info, save_path=None, file_path=None, timeout=60):
    """Return target HITRAN ``.par`` and ``.pf`` paths."""
    save_path = resolve_save_path(save_path=save_path, file_path=file_path)
    targets = []
    for molecule, isotopologue, _, _ in _normalize_hitran_configs(species_info, timeout=timeout):
        base_path = _target_base(save_path, molecule, isotopologue)
        targets.extend([base_path + '.par', base_path + '.pf'])
    return targets


def download_hitran(
    species_info,
    save_path=None,
    file_path=None,
    url_path=None,
    download=True,
    timeout=60,
):
    """Save HITRAN URLs and optionally download ``.par``/``.pf`` files."""
    save_path = resolve_save_path(save_path=save_path, file_path=file_path)
    if url_path is None:
        url_path = default_url_path(save_path, 'hitran')

    configs = _normalize_hitran_configs(species_info, timeout=timeout)
    urls = []
    for _, _, meta, wn_range in configs:
        urls.extend([_line_url(meta['global_id'], wn_range), meta['q_url']])
    _print_hitran_urls(configs)
    write_urls(urls, url_path)
    print('\nAll HITRAN URLs have been saved to', url_path)

    if not download:
        return urls

    for molecule, isotopologue, meta, wn_range in configs:
        table_name = f'{molecule}__{isotopologue}'
        output_dir = os.path.join(save_path, molecule, isotopologue)
        line_url = _line_url(meta['global_id'], wn_range)
        download_url_to_path(
            line_url,
            os.path.join(output_dir, table_name + '.par'),
            timeout=timeout,
            expected='text',
        )
        download_url_to_path(
            meta['q_url'],
            os.path.join(output_dir, table_name + '.pf'),
            timeout=timeout,
            expected='text',
        )
        print(f'HITRAN - {molecule} - {isotopologue}: downloaded .par and .pf')

    return urls
