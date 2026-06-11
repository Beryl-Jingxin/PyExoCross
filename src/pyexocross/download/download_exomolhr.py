"""Download ExoMolHR CSV line lists and partition functions."""

import os
import json
from urllib.parse import urlencode, urljoin, urlparse, parse_qsl, urlunparse

import requests
from bs4 import BeautifulSoup

from .common import (
    default_url_path,
    download_url_to_path,
    download_zip_csv_to_path,
    resolve_save_path,
    write_urls,
)


EXOMOLHR_MASTER_URL = 'https://www.exomol.com/exomolhr/exomolhr.all.json'
EXOMOLHR_DB_URL = 'https://www.exomol.com/exomolhr/db'
EXOMOLHR_GET_DATA_URL = 'https://www.exomol.com/exomolhr/get-data/'
EXOMOLHR_PF_URL = 'https://www.exomol.com/exomolhr/pf'


def _normalize_molecule_isotopologues(species_info):
    normalized = []
    for molecule, isotopologues in species_info.items():
        if isinstance(isotopologues, str):
            normalized.append((molecule, isotopologues, None))
            continue
        if isinstance(isotopologues, (list, tuple)):
            for isotopologue in isotopologues:
                normalized.append((molecule, isotopologue, None))
            continue
        if isinstance(isotopologues, dict):
            for isotopologue, config in isotopologues.items():
                normalized.append((molecule, isotopologue, config))
            continue
        raise ValueError(f'Invalid ExoMolHR isotopologue config for {molecule}: {isotopologues}')
    return normalized


def _online_csv_page_url(isotopologue, config):
    if not isinstance(config, dict):
        raise ValueError(f'ExoMolHR online query config for {isotopologue} must be a dict or None.')
    wn_range = config.get('wn_range')
    temperature = config.get('T', config.get('temperature'))
    threshold = config.get('threshold', config.get('intensity_threshold', config.get('Smin')))
    if wn_range is None or temperature is None or threshold is None:
        raise ValueError(
            f'ExoMolHR online query for {isotopologue} requires T, wn_range, and threshold.'
        )
    params = [
        ('iso', isotopologue),
        ('numin', wn_range[0]),
        ('numax', wn_range[1]),
        ('T', temperature),
        ('Smin', threshold),
    ]
    return EXOMOLHR_GET_DATA_URL + '?' + urlencode(params)


def _resolve_online_csv_download_url(page_url, timeout=60):
    response = requests.get(page_url, timeout=timeout)
    response.raise_for_status()

    content_type = response.headers.get('content-type', '').lower()
    content_disposition = response.headers.get('content-disposition', '').lower()
    if 'text/csv' in content_type or 'attachment' in content_disposition:
        return response.url

    soup = BeautifulSoup(response.text, 'html.parser')

    for link in soup.find_all('a', href=True):
        href = link.get('href')
        text = link.get_text(' ', strip=True).lower()
        href_lower = href.lower()
        if 'download' in text or 'download' in href_lower or href_lower.endswith('.csv'):
            return urljoin(response.url, href)

    for button in soup.find_all(['button', 'input']):
        text = button.get_text(' ', strip=True).lower()
        value = str(button.get('value', '')).lower()
        if 'download' not in text and 'download' not in value:
            continue
        form = button.find_parent('form')
        if form is None:
            continue
        action = form.get('action') or response.url
        form_url = urljoin(response.url, action)
        method = str(form.get('method', 'get')).lower()
        if method != 'get':
            return form_url
        params = []
        for field in form.find_all(['input', 'select', 'textarea']):
            name = field.get('name')
            if not name:
                continue
            value = field.get('value', '')
            params.append((name, value))
        if not params:
            return form_url
        parsed = urlparse(form_url)
        query = parse_qsl(parsed.query, keep_blank_values=True)
        query.extend(params)
        return urlunparse(parsed._replace(query=urlencode(query)))

    raise ValueError(
        'Could not resolve ExoMolHR online CSV download URL from page: '
        f'{page_url}'
    )


def _online_csv_filename(molecule, isotopologue, config):
    wn_range = config.get('wn_range')
    temperature = config.get('T', config.get('temperature'))
    threshold = config.get('threshold', config.get('intensity_threshold', config.get('Smin')))
    threshold_text = f'{threshold:g}' if isinstance(threshold, float) else str(threshold)
    threshold_text = threshold_text.replace('+', '')
    return (
        f'{molecule}__{isotopologue}'
        f'__T{temperature}'
        f'__{wn_range[0]}-{wn_range[1]}'
        f'__Smin{threshold_text}.csv'
    )


def _infer_dataset_from_master(master_data, molecule, isotopologue):
    molecule_data = master_data.get('molecules', {}).get(molecule)
    if molecule_data is None:
        raise ValueError(f'No ExoMolHR molecule found in master JSON: {molecule}')
    datasets = [
        row.get('dataset_name')
        for row in molecule_data.get('linelist', [])
        if row.get('iso_slug') == isotopologue
    ]
    datasets = [dataset for dataset in datasets if dataset]
    if len(set(datasets)) == 1:
        return datasets[0]
    raise ValueError(
        f'Cannot infer unique ExoMolHR dataset for {molecule}/{isotopologue}. '
        'The ExoMolHR API must expose a single dataset for each isotopologue.'
    )


def _get_master_data(timeout):
    response = requests.get(EXOMOLHR_MASTER_URL, timeout=timeout)
    response.raise_for_status()
    return json.loads(response.text)


def _exomolhr_entries(species_info, save_path=None, timeout=60, resolve_online=True):
    master_data = None
    entries = []
    for molecule, isotopologue, config in _normalize_molecule_isotopologues(species_info):
        if master_data is None:
            master_data = _get_master_data(timeout)
        dataset = _infer_dataset_from_master(master_data, molecule, isotopologue)
        if config is None:
            csv_url = f'{EXOMOLHR_DB_URL}/{molecule}__{isotopologue}__{dataset}.csv'
        else:
            page_url = _online_csv_page_url(isotopologue, config)
            csv_url = (
                _resolve_online_csv_download_url(page_url, timeout=timeout)
                if resolve_online
                else page_url
            )
        pf_url = f'{EXOMOLHR_PF_URL}/{isotopologue}__{dataset}.pf'
        if save_path is None:
            csv_path = None
            pf_path = None
        else:
            base_dir = os.path.join(save_path, molecule, isotopologue)
            if config is None:
                csv_filename = f'{molecule}__{isotopologue}.csv'
            else:
                csv_filename = _online_csv_filename(molecule, isotopologue, config)
            csv_path = os.path.join(base_dir, csv_filename)
            pf_path = os.path.join(base_dir, f'{isotopologue}__{dataset}.pf')
        entries.extend([(csv_url, csv_path), (pf_url, pf_path)])
    return entries


def get_exomolhr_urls(species_info, timeout=60, resolve_online=True):
    """Return ExoMolHR ``.csv`` and ``.pf`` URLs."""
    entries = _exomolhr_entries(
        species_info,
        timeout=timeout,
        resolve_online=resolve_online,
    )
    urls = [url for url, _ in entries]
    _print_exomolhr_entries(species_info, entries)
    return urls


def _print_exomolhr_entries(species_info, entries):
    """Print ExoMolHR URL groups in the same style as ExoMol downloads."""
    entry_index = 0
    for molecule, isotopologue, _ in _normalize_molecule_isotopologues(species_info):
        group_urls = [url for url, _ in entries[entry_index:entry_index + 2]]
        entry_index += 2
        print(f'ExoMolHR - {molecule} - {isotopologue}: {len(group_urls)} file(s)')
        for url in group_urls:
            print(url)


def download_exomolhr(
    species_info,
    save_path=None,
    file_path=None,
    url_path=None,
    download=True,
    timeout=60,
    resolve_online=True,
):
    """Save ExoMolHR URLs and optionally download them."""
    save_path = resolve_save_path(save_path=save_path, file_path=file_path)
    if url_path is None:
        url_path = default_url_path(save_path, 'exomolhr')

    entries = _exomolhr_entries(
        species_info,
        save_path=save_path,
        timeout=timeout,
        resolve_online=resolve_online,
    )
    urls = [url for url, _ in entries]
    _print_exomolhr_entries(species_info, entries)
    write_urls(urls, url_path)
    print('\nAll ExoMolHR URLs have been saved to', url_path)

    if download:
        for url, target_path in entries:
            if target_path.endswith('.csv') and '/get-data/download/' in url:
                download_zip_csv_to_path(url, target_path, timeout=timeout)
            else:
                expected = 'csv' if target_path.endswith('.csv') else 'text'
                download_url_to_path(url, target_path, timeout=timeout, expected=expected)
        print('\nAll ExoMolHR files have been downloaded to', save_path, 'folder!')

    return urls
