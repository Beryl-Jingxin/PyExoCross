# Import all what we need
import os
import json
import urllib3
import requests
import subprocess
from tqdm.notebook import tqdm
urllib3.disable_warnings()

# File Paths and Molecules
################## Could be changed ! ##################
url_path = '/home/jingxin/data/ExoMol/url/api__urls.txt'
file_path = '/home/jingxin/data/ExoMol/'
molecules = ['MgH']
########################################################

# Get API URLs
def get_api(molecules):
    molecule_str = []
    api_url = []
    for i in range(len(molecules)):
        molecule_str.append(molecules[i].replace('_p','+').split('__')[0].replace('+','_p'))
        api_url.append('http://exomol.com/api/?molecule=*&datatype=linelist'.replace('*',molecule_str[i]))
    return(api_url)

# Get Download Links with API
def get_urls(molecules):
    """Get the download url from API."""
    api_url = get_api(molecules)
    urls = []
    for i in tqdm(range(len(molecules))):
        response = requests.get(api_url[i])
        if(response.status_code != 200):
            print('ExoMol API Error' + str(response.status_code))

        # If the obtained status code is 200, it is correct.
        else:
            content = response.text            # Get the relevant content.
            json_dict = json.loads(content)    # Convert json into dictionary.
            iso_formulas = list(json_dict.keys())
            for iso_formula in iso_formulas:
                datasets = list(json_dict[iso_formula]['linelist'].keys())[1:]
                for dataset in datasets:
                    files_info = json_dict[iso_formula]['linelist'][dataset]
                    if files_info['recommended'] == True:
                        files_meta = files_info['files']
                        nfiles = len(files_meta)
                        trans_count = 0
                        trans_urls = []
                        for j in range(nfiles):
                            file_meta = files_meta[j]
                            url = "http://www." + file_meta.get('url')
                            if url.endswith('states.bz2'):
                                states_url = url
                                def_url = states_url.replace('.states.bz2','.def.json')
                                pf_url = states_url.replace('.states.bz2','.pf')
                            elif url.endswith('trans.bz2'):
                                trans_urls.append(url)
                                trans_count += 1
                            else:
                                print('No line list files.')
                        start = len(urls)
                        urls.extend([def_url, pf_url, states_url])
                        urls.extend(trans_urls)
                        print(f'{molecules[i]} - {iso_formula} - {dataset}: {trans_count} trans file(s)')
                        for entry in urls[start:]:
                            print(entry)
                        
    return(urls) 

# Download line list Files
# We write all the download URLs into a text file, name it as api__urls.txt. 
# In Linux, we use command:
# wget  -r -nH --cut-dirs=1 -P savePath -i PathOFapi__urls.txt
# Download line list files with urls and save them into correspoding folders.
def download_files(molecules, url_path):
    urls = get_urls(molecules)
    # Save all URLs to a text file
    if os.path.exists(url_path):
        pass
    else:
        os.makedirs(url_path, exist_ok=True)
    with open(url_path, "w", encoding="utf-8") as fh:
        for entry in urls:
            fh.write(f"{entry}\n")
    print('\nAll URLs have been saved to', url_path+'api__urls.txt')
    command = f'wget -r -nH --cut-dirs=1 -P {file_path} -i {url_path}'
    subprocess.run(command, shell=True)
    print('\nAll files have been downloaded to', file_path, 'folder!')

download_files(molecules, url_path)                 