"""Download ExoAtom files from the ExoMol ExoAtom service."""

import os

from .common import default_url_path, download_url_file, resolve_save_path, write_urls


EXOATOM_BASE_URL = 'https://www.exomol.com/exoatom/db'
EXOATOM_EXTENSIONS = ['.adef.json', '.states', '.trans', '.pf']


def _dataset_from_config(config):
    if isinstance(config, dict):
        return config.get('dataset')
    return config


def _normalize_atom_datasets(species_info):
    normalized = []
    for atom, config in species_info.items():
        if isinstance(config, str):
            normalized.append((atom, None, config))
            continue
        if not isinstance(config, dict):
            raise ValueError(f'Invalid ExoAtom config for {atom}: {config}')

        dataset = config.get('dataset')
        if dataset is not None:
            normalized.append((atom, None, dataset))

        for isotopologue, isotope_config in config.items():
            if isotopologue == 'dataset':
                continue
            isotope_dataset = _dataset_from_config(isotope_config)
            if isotope_dataset is None:
                raise ValueError(f'Missing dataset for ExoAtom {atom} isotopologue {isotopologue}.')
            normalized.append((atom, isotopologue, isotope_dataset))
    return normalized


def get_exoatom_urls(species_info):
    """Return ExoAtom URLs for selected atoms/isotopologues/datasets."""
    urls = []
    for atom, isotopologue, dataset in _normalize_atom_datasets(species_info):
        if isotopologue is None:
            base_url = f'{EXOATOM_BASE_URL}/{atom}/{dataset}/{atom}__{dataset}'
            label = atom
        else:
            base_url = f'{EXOATOM_BASE_URL}/{atom}/{isotopologue}/{dataset}/{isotopologue}__{dataset}'
            label = isotopologue
        group_urls = [base_url + extension for extension in EXOATOM_EXTENSIONS]
        urls.extend(group_urls)
        print(f'ExoAtom - {label} - {dataset}: {len(group_urls)} file(s)')
        for url in group_urls:
            print(url)
    return urls


def download_exoatom(
    species_info,
    save_path=None,
    file_path=None,
    url_path=None,
    download=True,
):
    """Save ExoAtom URLs and optionally download them."""
    save_path = resolve_save_path(save_path=save_path, file_path=file_path)
    if url_path is None:
        url_path = default_url_path(save_path, 'exoatom')

    urls = get_exoatom_urls(species_info)
    write_urls(urls, url_path)
    print('\nAll ExoAtom URLs have been saved to', url_path)

    if download:
        download_url_file(url_path, save_path)
        print('\nAll ExoAtom files have been downloaded to', save_path, 'folder!')

    return urls
