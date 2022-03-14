# Import all what we need
import os
import sys
import json
import requests
import pandas as pd
from io import StringIO
from tqdm import tqdm


# Change the Information Here
abspath = '/home/jingxin/'


# Read Master File
# Read the master file 'ExoMol All'. The URL is http://www.exomol.com/db/exomol.all.
exomol_all_url = 'http://www.exomol.com/db/exomol.all'
content = requests.get(exomol_all_url).text.replace('#','')
exomol_col_name = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6']
exomol_all = pd.read_csv(StringIO(content), sep='\\s+', names=exomol_col_name, header=None)


# Get all molecules, their iso-slugs, isoFormula and isotopologue dataset names. 
first = exomol_all['c1']
fourth = exomol_all['c4']
row = len(first)
iso_slug = pd.DataFrame()
iso_formula = pd.DataFrame()
iso_dataset = pd.DataFrame()
molecule_unduplicated = pd.DataFrame()
num_iso_in_mol = pd.DataFrame()

for i in tqdm(range(row)):
    _iso_slug = exomol_all[first.isin(['Iso-slug'])]['c0'].values
    _iso_formula = exomol_all[first.isin(['IsoFormula'])]['c0'].values
    _iso_dataset = exomol_all[first.isin(['Isotopologue'])]['c0'].values
    _molecule_unduplicated = exomol_all[first.isin(['Molecule'])]['c0'].values
    _num_iso_in_mol = exomol_all[fourth.isin(['considered'])]['c0'].values

iso_slug = iso_slug.append(pd.DataFrame(_iso_slug))
iso_formula = iso_formula.append(pd.DataFrame(_iso_formula))
iso_dataset = iso_dataset.append(pd.DataFrame(_iso_dataset))
molecule_unduplicated = molecule_unduplicated.append(pd.DataFrame(_molecule_unduplicated))
num_iso_in_mol = num_iso_in_mol.append(pd.DataFrame(_num_iso_in_mol))
isotopologue_unduplicated = iso_slug.drop_duplicates(subset=None, keep='first', inplace=False)
iso_datset_unduplicated = iso_dataset.drop_duplicates(subset=None, keep='first', inplace=False)


# Set the length of molecules list to be the same as the length of iso-slugs and isotopologue dataset names for following loop.
molecule_repeated = pd.DataFrame()
molecule_num = len(molecule_unduplicated)

for j in range(molecule_num):    
    molecule_repeated = molecule_repeated.append(pd.DataFrame((molecule_unduplicated.values[j] + ' ')
                                                              * int(num_iso_in_mol.values[j])))

molecule_str = (str(molecule_repeated.values).replace("[['"," ")
                .replace("']\n ['"," ").replace("']]"," ").replace("+","_p"))
molecule = pd.read_csv(StringIO(molecule_str), sep='\s+', header=None)


# Read Def File
# Get all URLs of def files. The number of def files should be the same as the number of isotopologue datasets. 
# The URLs contains the names of molecules, iso-slugs and isotopologue datasets. 
# We save their corresponding isoFormula names as another column.
def_url = pd.DataFrame()
def_num = len(iso_slug)
for i in tqdm(range(def_num)):
    def_url = def_url.append('http://www.exomol.com/db/' + molecule[i] + '/'
                             + iso_slug.values[i] + '/'+ iso_dataset.values[i] + '/'
                             + iso_slug.values[i] + '__' + iso_dataset.values[i] + '.def')

def_url['IsoFormula'] = iso_formula
def_url.columns = ['def url','IsoFormula']


# Download def files and save them into ./data/def/ folder. 
# Save the names of these def files with all information we got before, that is to say, save as 'molecule_isoFormula_iso-slug_isotopologue.def'. 
# It will be more convenient for processing data later.
def download_deffile(file_url):
    
    failed_list = [] 
    for _link in tqdm(file_url['def url'].values):
        
        link = _link
        iso_slug = link.split('/')[-4]
        iso_formula = (str(file_url[file_url['def url'].isin([link])]['IsoFormula'].values)
                       .replace("['","").replace("']",""))
        inital_def_name = link.split('/')[-1]
        new_def_filename = iso_slug + '_' + iso_formula + '_' + inital_def_name
        print("Downloading file: %s" % new_def_filename)
        print(link)
 
        # Make folders for save doanloaded files.
        folder_name = abspath + '/data/def/'
        if os.path.exists(folder_name):
            pass
        else:
            os.makedirs(folder_name, exist_ok=True)
        filename = os.path.join(folder_name, new_def_filename)
        
        try:
            r = requests.get(link, stream=True, verify=False)
        except Exception:
            failed_list.append(new_def_filename)
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
    return

download_deffile(def_url)


# Extract the def files whose the uncertainty availability row shows YES.
path = abspath + '/data/def'
def_col_name = ['c0', '#', 'c1', 'c2', 'c3', 'c4', 'c5']

for(dirpath,dirnames,files)in os.walk(path):
    tot = 0
    count = 0
    unc_def_filename = []
    for filename in files:
        filepath = os.path.join(dirpath, filename)
        tot += 1
        def_df =  pd.read_csv(filepath,sep='\s+', names=def_col_name, header=None)
        c1 = def_df['c1']
        if def_df[c1.isin(['Uncertainty'])]['c0'].values == '1':
            unc_def_filename.append(filename)
            count += 1            
        else:
            no_unc_def_filename = []
        
    print('There are ', tot, ' def files.\n')
    print('The uncertainty availability shows YES in the following ', count, ' def files:\n', unc_def_filename)
    print('\nThe uncertainty availability does not exit or shows NO in other ', tot - count, 'def files.')


# Get Download Links with API
# Get the API URLs of those uncertainty molecules.
molecule_str = []
iso_formula_str = []
iso_dataset_str = []
iso = []
api_url = []
target_link = []

for i in range(len(unc_def_filename)):
    molecule_str.append(unc_def_filename[i].replace('_p','+').split('_')[0].replace('+','_p'))
    iso_formula_str.append(unc_def_filename[i].replace('_p','+').split('_')[1])
    iso_dataset_str.append(unc_def_filename[i].split('_')[-1].split('.')[0])
    
    _iso = (iso_formula_str[i], iso_dataset_str[i])
    iso.append(_iso)
    
    api_url.append('http://exomol.com/api/?molecule=*&datatype=linelist'.replace('*',molecule_str[i]))
    
print(iso)


# Get the download links of states.bz2 files and trans.bz2 files from API.
def get_target_url():
    """
    Get the download url from API.

    """
    file_url = []
    for i in range(len(iso)):
        response = requests.get(api_url[i])
        if(response.status_code != 200):
            print('ExoMol API Error' + str(response.status_code))

        # If the obtained status code is 200, it is correct.
        else:
            content = response.text            # Get the relevant content.
            json_dict = json.loads(content)    # Convert json into dictionary.

            # Extract files information from dictionary and convert them into list
            iso_slug = iso[i][0]
            iso_dataset = iso[i][1]
            json_list = json_dict[iso_slug]['linelist'][iso_dataset]['files']
            print('The number of downloading files for', iso_slug, iso_dataset, ': ', len(json_list))
            print("Download links:")

            url_show = []
            for j in range(len(json_list)):
                link = json_list[j].get('url')
                try:
                    if((link.endswith('states.bz2') or link.endswith('trans.bz2'))):
                        file_url.append("http://www." + link)
                        url_show.append("http://www." + link)
                except KeyError:
                    print('Keyerror, keep going!')
                    
        for k in url_show:
            print(k)

    return file_url

target_link = get_target_url()


# Download States and Trans Files
# We write all the download URLs into a text file, name it as $api__urls.txt$. 
# In Linux, we use command 'wget -d -r -i /.../save_path/.../api__urls.txt'
# Download states.bz2 files and trans.bz2 files with download links. Save these files into correspoding folders.
url_path = abspath + '/data/url'
url_filename = url_path + '/api__urls.txt'

if os.path.exists(url_path):
    pass
else:
    os.makedirs(url_path, exist_ok=True)

with open(url_filename, 'w') as file:
    file.write('\n'.join(target_link))
