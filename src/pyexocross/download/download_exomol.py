"""Download recommended ExoMol line-list files through the ExoMol API."""

import os
import re
import json
import requests
from tqdm import tqdm

from .common import default_url_path, download_url_file, resolve_save_path, write_urls


def _api_url(molecule):
    molecule_name = molecule.replace('_p', '+').split('__')[0].replace('+', '_p')
    return f'https://exomol.com/api/?molecule={molecule_name}&datatype=linelist'


def _normalize_molecule_isotopologues(species_info):
    molecules = list(species_info.keys())
    isotopologue_configs = []
    for molecule in molecules:
        molecule_isos = species_info[molecule]
        if molecule_isos in (None, ''):
            isotopologue_configs.append(None)
        elif isinstance(molecule_isos, str):
            isotopologue_configs.append({molecule_isos: {'wn_range': None}})
        elif isinstance(molecule_isos, dict):
            isotopologue_config = {}
            for isotopologue, config in molecule_isos.items():
                if isinstance(config, dict):
                    isotopologue_config[isotopologue] = config
                else:
                    isotopologue_config[isotopologue] = {'wn_range': config}
            isotopologue_configs.append(isotopologue_config)
        else:
            isotopologue_configs.append(
                {isotopologue: {'wn_range': None} for isotopologue in molecule_isos}
            )
    return molecules, isotopologue_configs


def _get_wn_range(isotopologue_config):
    if isotopologue_config is None:
        return None
    return isotopologue_config.get('wn_range')


def _strict_states_filename(isotopologue, dataset):
    return f'{isotopologue}__{dataset}.states.bz2'


def _strict_trans_filename(isotopologue, dataset):
    return f'{isotopologue}__{dataset}.trans.bz2'


def _strict_segmented_trans_pattern(isotopologue, dataset):
    return re.compile(
        rf'^{re.escape(isotopologue)}__{re.escape(dataset)}__(\d+)-(\d+)\.trans\.bz2$'
    )


def _infer_iso_slug_from_url(url, dataset):
    filename = os.path.basename(url)
    suffixes = [
        f'__{dataset}.states.bz2',
        f'__{dataset}.trans.bz2',
    ]
    for suffix in suffixes:
        if filename.endswith(suffix):
            return filename[:-len(suffix)]

    segmented_match = re.match(rf'^(.+)__{re.escape(dataset)}__\d+-\d+\.trans\.bz2$', filename)
    if segmented_match is not None:
        return segmented_match.group(1)

    return None


def _trans_url_in_wn_range(url, isotopologue, dataset, wn_range):
    filename = os.path.basename(url)
    if filename == _strict_trans_filename(isotopologue, dataset):
        return True

    match = _strict_segmented_trans_pattern(isotopologue, dataset).match(filename)
    if match is None:
        return False

    if wn_range in (None, []):
        return True

    wn_min, wn_max = wn_range
    file_min = int(match.group(1))
    file_max = int(match.group(2))
    return file_min >= wn_min and file_max <= wn_max


def get_exomol_urls(species_info, timeout=60, show_progress=True):
    """
    Get recommended ExoMol file URLs for selected molecules/isotopologues.

    Parameters
    ----------
    species_info : dict
        Mapping of molecule names to isotopologue configuration, for example::

            {
                'MgH': {
                    '24Mg-1H': {'wn_range': None},
                },
                'H2O': {
                    '1H2-16O': {'wn_range': [41000, 41200]},
                },
            }

        For segmented ``.trans.bz2`` files, ``wn_range`` keeps only files fully
        inside the requested range. Exact ``iso__dataset.trans.bz2`` files are
        kept because they cannot be range-filtered by filename.
    molecule_isotopologues : dict, optional
        Alias for ``species_info`` (backward compatibility).
    timeout : int or float, optional
        Timeout in seconds for each ExoMol API request.
    show_progress : bool, optional
        Show a tqdm progress bar.

    Returns
    -------
    list of str
        URLs for ``.def.json``, ``.pf``, ``.states.bz2``, and selected
        ``.trans.bz2`` files.
    """
    molecules, isotopologue_configs = _normalize_molecule_isotopologues(species_info)
    iterator = range(len(molecules))
    if show_progress:
        iterator = tqdm(iterator)

    urls = []
    for i in iterator:
        molecule = molecules[i]
        target_isotopologue_config = isotopologue_configs[i]
        response = requests.get(_api_url(molecule), timeout=timeout)
        if response.status_code != 200:
            print('ExoMol API Error' + str(response.status_code))
            continue

        json_dict = json.loads(response.text)
        found_isotopologues = set()
        for iso_formula, iso_info in json_dict.items():
            linelist_info = iso_info.get('linelist', {})
            for dataset, files_info in linelist_info.items():
                if not isinstance(files_info, dict):
                    continue
                if files_info.get('recommended') is not True:
                    continue

                trans_count = 0
                trans_urls = []
                states_url = None
                iso_slug = None
                for file_meta in files_info.get('files', []):
                    url = 'https://www.' + file_meta.get('url')
                    filename = os.path.basename(url)
                    inferred_iso_slug = _infer_iso_slug_from_url(url, dataset)
                    if inferred_iso_slug is not None:
                        iso_slug = inferred_iso_slug
                    if iso_slug is None:
                        continue
                    if target_isotopologue_config is not None and iso_slug not in target_isotopologue_config:
                        continue

                    current_wn_range = (
                        None if target_isotopologue_config is None
                        else _get_wn_range(target_isotopologue_config[iso_slug])
                    )
                    if filename == _strict_states_filename(iso_slug, dataset):
                        states_url = url
                        def_url = states_url.replace('.states.bz2', '.def.json')
                        pf_url = states_url.replace('.states.bz2', '.pf')
                    elif _trans_url_in_wn_range(url, iso_slug, dataset, current_wn_range):
                        trans_urls.append(url)
                        trans_count += 1

                if states_url is None:
                    if target_isotopologue_config is None or iso_slug in target_isotopologue_config:
                        print(f'{molecule} - {iso_slug or iso_formula} - {dataset}: no strict states file found.')
                    continue

                found_isotopologues.add(iso_slug)
                start = len(urls)
                urls.extend([def_url, pf_url, states_url])
                urls.extend(trans_urls)
                print(f'{molecule} - {iso_slug} - {dataset}: {trans_count} trans file(s)')
                for entry in urls[start:]:
                    print(entry)

        if target_isotopologue_config is not None:
            missing_isotopologues = [
                isotopologue for isotopologue in target_isotopologue_config
                if isotopologue not in found_isotopologues
            ]
            for missing_isotopologue in missing_isotopologues:
                print(f'{molecule} - {missing_isotopologue}: recommended line list not found in ExoMol API response.')

    return urls


def download_exomol(
    species_info,
    save_path=None,
    file_path=None,
    url_path=None,
    download=True,
    timeout=60,
    show_progress=True,
):
    """
    Save ExoMol API URLs and optionally download the selected line-list files.

    Parameters
    ----------
    species_info : dict
        Molecule/isotopologue configuration in dict mode.
    save_path : str
        Folder where downloaded files are saved when ``download=True``.
    file_path : str, optional
        Alias for ``save_path``.
    url_path : str, optional
        Path to the generated URL list. Defaults to
        ``<save_path>/url/api__urls.txt``.
    download : bool, optional
        If True, run ``wget`` with the generated URL list. If False, only save
        and return URLs.
    timeout : int or float, optional
        Timeout in seconds for each ExoMol API request.
    show_progress : bool, optional
        Show a tqdm progress bar while querying the API.

    Returns
    -------
    list of str
        The generated URL list.
    """
    save_path = resolve_save_path(save_path=save_path, file_path=file_path)
    if url_path is None:
        url_path = default_url_path(save_path, 'exomol')

    urls = get_exomol_urls(
        species_info,
        timeout=timeout,
        show_progress=show_progress,
    )
    write_urls(urls, url_path)
    print('\nAll URLs have been saved to', url_path)

    if download:
        download_url_file(url_path, save_path)
        print('\nAll files have been downloaded to', save_path, 'folder!')

    return urls
