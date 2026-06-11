# Import all what we need
import os
import re
import json
import urllib3
import requests
import subprocess
from tqdm import tqdm
urllib3.disable_warnings()

# File Paths and Molecules
################## Could be changed ! ##################
# Directory that will hold the generated api__urls.txt file
url_dir = '/Users/beryl/Academic/UCL/PhD/Data/database/ExoMol/url/'
# Full path to the urls file (derived from url_dir)
url_path = os.path.join(url_dir, 'api__urls.txt')
file_path = '/Users/beryl/Academic/UCL/PhD/Data/database/ExoMol/'
molecule_isotopologues = {
    'MgH': {
        '24Mg-1H': {'wn_range': None},
        '25Mg-1H': {'wn_range': None},
    },
    'NO': {
        '14N-16O': {'wn_range': None},
    },
    'H2O': {
        '1H2-16O': {'wn_range': [41000, 41200]},
    },
}
########################################################

# Get API URLs
def get_api(molecules):
    molecule_str = []
    api_url = []
    for i in range(len(molecules)):
        molecule_str.append(molecules[i].replace('_p','+').split('__')[0].replace('+','_p'))
        api_url.append('https://exomol.com/api/?molecule=*&datatype=linelist'.replace('*',molecule_str[i]))
    return(api_url)


# Get Download Links with API
def normalize_molecule_isotopologues(molecule_isotopologues):
    molecules = list(molecule_isotopologues.keys())
    isotopologue_configs = []
    for molecule in molecules:
        molecule_isos = molecule_isotopologues[molecule]
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


def get_wn_range(isotopologue_config):
    if isotopologue_config is None:
        return None
    return isotopologue_config.get('wn_range')


def strict_states_filename(isotopologue, dataset):
    return f'{isotopologue}__{dataset}.states.bz2'


def strict_trans_filename(isotopologue, dataset):
    return f'{isotopologue}__{dataset}.trans.bz2'


def strict_segmented_trans_pattern(isotopologue, dataset):
    return re.compile(
        rf'^{re.escape(isotopologue)}__{re.escape(dataset)}__(\d+)-(\d+)\.trans\.bz2$'
    )


def infer_iso_slug_from_url(url, dataset):
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


def trans_url_in_wn_range(url, isotopologue, dataset, wn_range):
    filename = os.path.basename(url)
    if filename == strict_trans_filename(isotopologue, dataset):
        return True

    match = strict_segmented_trans_pattern(isotopologue, dataset).match(filename)
    if match is None:
        return False

    if wn_range in (None, []):
        return True

    wn_min, wn_max = wn_range
    file_min = int(match.group(1))
    file_max = int(match.group(2))
    return file_min >= wn_min and file_max <= wn_max


def get_urls(molecule_isotopologues):
    """Get the download url from API."""
    molecules, isotopologue_configs = normalize_molecule_isotopologues(molecule_isotopologues)
    api_url = get_api(molecules)
    urls = []
    for i in tqdm(range(len(molecules))):
        target_isotopologue_config = isotopologue_configs[i]
        response = requests.get(api_url[i], timeout=60)
        if(response.status_code != 200):
            print('ExoMol API Error' + str(response.status_code))

        # If the obtained status code is 200, it is correct.
        else:
            content = response.text            # Get the relevant content.
            json_dict = json.loads(content)    # Convert json into dictionary.
            found_isotopologues = set()
            for iso_formula, iso_info in json_dict.items():
                linelist_info = iso_info.get('linelist', {})
                for dataset, files_info in linelist_info.items():
                    if not isinstance(files_info, dict):
                        continue
                    if files_info.get('recommended') == True:
                        files_meta = files_info.get('files', [])
                        nfiles = len(files_meta)
                        trans_count = 0
                        trans_urls = []
                        states_url = None
                        iso_slug = None
                        current_wn_range = None
                        for j in range(nfiles):
                            file_meta = files_meta[j]
                            url = "https://www." + file_meta.get('url')
                            filename = os.path.basename(url)
                            inferred_iso_slug = infer_iso_slug_from_url(url, dataset)
                            if inferred_iso_slug is not None:
                                iso_slug = inferred_iso_slug
                            if iso_slug is None:
                                continue
                            if target_isotopologue_config is not None and iso_slug not in target_isotopologue_config:
                                continue
                            current_wn_range = (
                                None if target_isotopologue_config is None
                                else get_wn_range(target_isotopologue_config[iso_slug])
                            )
                            if filename == strict_states_filename(iso_slug, dataset):
                                states_url = url
                                def_url = states_url.replace('.states.bz2','.def.json')
                                pf_url = states_url.replace('.states.bz2','.pf')
                            elif trans_url_in_wn_range(url, iso_slug, dataset, current_wn_range):
                                trans_urls.append(url)
                                trans_count += 1
                        if states_url is None:
                            if target_isotopologue_config is None or iso_slug in target_isotopologue_config:
                                print(f'{molecules[i]} - {iso_slug or iso_formula} - {dataset}: no strict states file found.')
                            continue
                        found_isotopologues.add(iso_slug)
                        start = len(urls)
                        urls.extend([def_url, pf_url, states_url])
                        urls.extend(trans_urls)
                        print(f'{molecules[i]} - {iso_slug} - {dataset}: {trans_count} trans file(s)')
                        for entry in urls[start:]:
                            print(entry)
            if target_isotopologue_config is not None:
                missing_isotopologues = [
                    iso for iso in target_isotopologue_config
                    if iso not in found_isotopologues
                ]
                for missing_isotopologue in missing_isotopologues:
                    print(f'{molecules[i]} - {missing_isotopologue}: recommended line list not found in ExoMol API response.')
                        
    return(urls) 

# Download line list Files
# We write all the download URLs into a text file, name it as api__urls.txt. 
# In Linux, we use command:
# wget  -r -nH --cut-dirs=1 -P savePath -i PathOFapi__urls.txt
# Download line list files with urls and save them into correspoding folders.
def download_files(molecule_isotopologues, url_path):
    urls = get_urls(molecule_isotopologues)
    # Save all URLs to a text file
    os.makedirs(os.path.dirname(url_path), exist_ok=True)
    with open(url_path, "w", encoding="utf-8") as fh:
        for entry in urls:
            fh.write(f"{entry}\n")
    print('\nAll URLs have been saved to', url_path)
    command = f'wget -r -nH --cut-dirs=1 -P {file_path} -i {url_path}'
    subprocess.run(command, shell=True)
    print('\nAll files have been downloaded to', file_path, 'folder!')

if __name__ == '__main__':
    download_files(molecule_isotopologues, url_path)
