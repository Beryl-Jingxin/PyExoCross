# Import all what we need
import os
import sys
import json
import urllib3
import requests
import subprocess
from tqdm import notebook
urllib3.disable_warnings()

# File Paths
################## Could be changed ! ##################
def_path = '/home/jingxin/data/test/def'
url_path = '/home/jingxin/data/pyexocross/url'
molecule = ['MgH','NO']
iso_slug = ['24Mg-1H','14N-16O']
iso_formula = ['(24Mg)(1H)','(14N)(16O)']
dataset = ['XAB','XABC']
########################################################

# Read Def File
# Get all URLs of def files. The number of def files should be the same as the number of isotopologue datasets. 
# The URLs contains the names of molecules, iso-slugs and isotopologue datasets. 
# We save their corresponding isoFormula names as another column.
def read_deffile(molecule, iso_slug, dataset):
    def_url = []
    def_num = len(iso_slug)
    for i in notebook.tqdm(range(def_num)):
        url = ('http://www.exomol.com/db/' + molecule[i] + '/'
               + iso_slug[i] + '/'+ dataset[i] + '/'
               + iso_slug[i] + '__' + dataset[i] + '.def')
        def_url.append(url)
    return(def_url)

# Download def files and save them into ./data/def/ folder. 
# Save the names of these def files with all information we got before, 
# that is to say, save as 'molecule_isoFormula_iso-slug_isotopologue.def'. 
# It will be more convenient for processing data later.
def download_deffile(def_path):
    failed_list = [] 
    def_url = read_deffile(molecule, iso_slug, dataset)
    for link in tqdm(def_url):
        def_filename = link.split('/')[-1]
        print("Downloading file: %s" % def_filename)
        print(link)
 
        # Make folders for save doanloaded files.
        if os.path.exists(def_path):
            pass
        else:
            os.makedirs(def_path, exist_ok=True)
        filename = os.path.join(def_path, def_filename)
        
        try:
            r = requests.get(link, stream=True, verify=False)
        except Exception:
            failed_list.append(def_filename)
            print(' download failed. Go to download next one\n')
              
        # For compute the progess.
        total_size = int(r.headers['Content-Length'])
        temp_size = 0    
   
        # Download started.
        with open(filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    temp_size += len(chunk)
                    f.write(chunk)
                    f.flush()
                    done = int(50 * temp_size / total_size)
                    sys.stdout.write("\r[%s%s] %d%%" % ('█' * done, ' ' * (50 - done),
                                                        100 * temp_size / total_size))
                    sys.stdout.flush()
        print(" Downloaded!\n")
    print("All def files downloaded!\n")    
    print("The files which are failed to download: \n")
    print(failed_list) # Record which file is failed to download.

# Get Download Links with API
# Get the API URLs of those uncertainty molecules.
def get_api(def_path):
    molecule_str = []
    iso_str = []
    api_url = []
    for i in range(len(molecule)):
        molecule_str.append(molecule[i].replace('_p','+').split('__')[0].replace('+','_p'))
        iso_str.append(iso_slug[i].replace('_p','+'))
        api_url.append('http://exomol.com/api/?molecule=*&datatype=linelist'.replace('*',molecule_str[i]))
    return(api_url)

# Get the download links of states.bz2 files and trans.bz2 files from API.
def get_target_url(def_path):
    """Get the download url from API."""
    file_url = []
    api_url = get_api(def_path)
    for i in range(len(molecule)):
        response = requests.get(api_url[i])
        if(response.status_code != 200):
            print('ExoMol API Error' + str(response.status_code))

        # If the obtained status code is 200, it is correct.
        else:
            content = response.text            # Get the relevant content.
            json_dict = json.loads(content)    # Convert json into dictionary.

            # Extract files information from dictionary and convert them into list
            _iso_formula = iso_formula[i]
            _dataset = dataset[i]
            json_list = json_dict[_iso_formula]['linelist'][_dataset]['files']
            onefile = "http://www." + json_list[0].get('url')
            def_pf = [onefile.replace('.states.bz2','.def'), onefile.replace('.states.bz2','.pf')]
            url_show = []
            for j in range(len(json_list)):
                link = json_list[j].get('url')
                try:
                    if((link.endswith('states.bz2') or link.endswith('trans.bz2'))):
                        print(link.split('/')[-1].replace('.states.bz2','.def'))
                        file_url.append("http://www." + link)
                        url_show.append("http://www." + link)
                except KeyError:
                    print('Keyerror, keep going!')
            file_url += def_pf
            url_show += def_pf
        print('\nThe number of downloading files for', _iso_formula, _dataset, ': ', len(url_show))
        print("Download links:")                    
        for k in url_show:
            print(k)
    return (file_url)

# Download States and Trans Files
# We write all the download URLs into a text file, name it as api__urls.txt. 
# In Linux, we use command:
# wget -d -r -i /.../save_path/.../api__urls.txt
# Download states.bz2 files and trans.bz2 files with download links. Save these files into correspoding folders.
def download_files(url_path):
    url_filename = url_path + '/api__urls.txt'

    if os.path.exists(url_path):
        pass
    else:
        os.makedirs(url_path, exist_ok=True)

    target_link = get_target_url(def_path)
    with open(url_filename, 'w') as file:
        file.write('\n'.join(target_link))
    command = f'wget -d -r -i {url_filename}'
    subprocess.run(command, shell=True)
    
download_deffile(def_path)
download_files(url_path)
