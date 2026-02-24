# Examples

This page contains complete, ready-to-run examples for every supported
database and function.  Replace the paths with your own data locations.

---

## ExoMol Examples

### Conversion (ExoMol -> HITRAN)

```python
import pyexocross as px

px.conversion(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    species_id=501,
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/', 
    logs_path='/path/to/output/log/exomol_conversion.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['par', 'e/f', 'eS', 'v', 'Lambda', 'Sigma', 'Omega'],
    qnsformat_list=['%1s', '%1s', '%13s', '%3d', '%2d', '%7.1f', '%7.1f'],
    conversion_format=1,
    conversion_min_freq=0,
    conversion_max_freq=30000,
    conversion_unc=0.01, 
    conversion_threshold=1e-30,
    global_qn_label_list=['eS', 'v', 'Omega'],
    global_qn_format_list=['%9s', '%2d', '%4s'],
    local_qn_label_list=['J', 'e/f'],
    local_qn_format_list=['%5.1f', '%2s']
)
```

### Partition Functions

```python
import pyexocross as px

px.partition_functions(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_pf.log',
    ntemp=1,
    tmax=5000
)
```

### Specific Heats

```python
import pyexocross as px

px.specific_heats(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_cp.log',
    ntemp=1,
    tmax=5000
)
```

### Cooling Functions

```python
import pyexocross as px

px.cooling_functions(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_cf.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    ntemp=1,
    tmax=5000
)
```

### Lifetimes

```python
import pyexocross as px

px.lifetimes(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_lifetime.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    compress=False
)
```

### Oscillator Strengths

```python
import pyexocross as px

px.oscillator_strengths(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_os.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    gf_or_f='f',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN',
    plot_unit='cm-1',
    limit_yaxis=1e-30
)
```

### Stick Spectra (LTE)

```python
import pyexocross as px

px.stick_spectra(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_stick.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    temperatures=[1000, 2000],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    abs_emi='Ab',
    plot=True,
    plot_method='log'
)
```

### Stick Spectra (Non-LTE, Treanor Method)

```python
import pyexocross as px

px.stick_spectra(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_stick_nlte.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['par', 'e/f', 'eS', 'v', 'Lambda', 'Sigma', 'Omega'],
    qnsformat_list=['%1s', '%1s', '%13s', '%3d', '%2d', '%7.1f', '%7.1f'],
    nlte_method='T',
    tvib_list=[1000, 2000, 3000],
    trot_list=[100, 200],
    vib_label=['v', 'eS'],
    rot_label=['J', 'e/f'],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    abs_emi='Ab',
    unc_filter=0.01,
    threshold=1e-30,
    qns_filter={
        'par': [],
        'e/f': [],
        'eS': [],
        'v': ['0,', '1,', '2,', '3,', '4,', ',0', ',1', ',2', ',3', ',4'],
    }
)
```

### Cross Sections

```python
import pyexocross as px

px.cross_sections(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/Databases/ExoMol/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol_xsec.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    temperatures=[1000, 2000, 3000],
    pressures=[1.0, 5.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    bin_size=0.1,
    profile='SciPyVoigt',
    predissociation=False,
    broadeners=['Default'],
    ratios=[1.0],
    cutoff=25.0,
    abs_emi='Ab',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN',
    plot_unit='cm-1',
    limit_yaxis=1e-30
)
```

### Cross Sections (Non-LTE, Treanor Method)

```python
import pyexocross as px

px.cross_sections(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    # read_path='/path/to/Databases/ExoMol/',
    # save_path='/path/to/output/',
    # logs_path='/path/to/output/log/exomol_xsec_nlte.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['par', 'e/f', 'eS', 'v', 'Lambda', 'Sigma', 'Omega'],
    qnsformat_list=['%1s', '%1s', '%13s', '%3d', '%2d', '%7.1f', '%7.1f'],
    nlte_method='T',
    tvib_list=[1000, 2000, 3000],
    trot_list=[100, 200],
    vib_label=['v', 'eS'],
    rot_label=['J', 'e/f'],
    # temperatures=[1000, 2000, 3000],
    pressures=[1.0, 5.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    bin_size=0.1,
    profile='SciPyVoigt',
    predissociation=False,
    broadeners=['Default'],
    ratios=[1.0],
    cutoff=25.0,
    abs_emi='Em',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN',
    plot_unit='cm-1',
    limit_yaxis=1e-30,
    read_path='/Users/beryl/Academic/UCL/PhD/Data/database/ExoMol/',
    save_path='/Users/beryl/Academic/UCL/PhD/Data/pyexocross/',
    logs_path='/Users/beryl/Academic/UCL/PhD/Data/pyexocross/log/test_api_exomol.log',
)
```

---

## ExoAtom Examples

### Conversion (ExoAtom -> HITRAN)

```python
import pyexocross as px

px.conversion(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exomol.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['configuration', 'Multiple', 'parity'],
    qnsformat_list=['%20s', '%10s', '%2s'],
    conversion_format=1,
    conversion_min_freq=0,
    conversion_max_freq=115400,
    conversion_unc=None, 
    conversion_threshold=None,
    global_qn_label_list=['configuration'],
    global_qn_format_list=['%20s'],
    local_qn_label_list=['J', 'Multiple', 'parity'],
    local_qn_format_list=['%5.1f', '%10s', '%1s'],
)
```

### Partition Functions

```python
import pyexocross as px

px.partition_functions(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exoatom_pf.log',
    ntemp=1,
    tmax=6000,
)
```

### Specific Heats

```python
import pyexocross as px

px.specific_heats(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exoatom_cp.log',
    ntemp=1,
    tmax=6000,
)
```

### Cooling Functions

```python
import pyexocross as px

px.cooling_functions(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exoatom_cf.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    ntemp=1,
    tmax=6000,
)
```

### Lifetimes

```python
import pyexocross as px

px.lifetimes(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exoatom_lifetime.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    compress=False,
)
```

### Oscillator Strengths

```python
import pyexocross as px

px.oscillator_strengths(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exoatom_os.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    gf_or_f='f',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN',   
    plot_unit='cm-1',     
    limit_yaxis=1e-30, 
)
```

### Stick Spectra (Non-LTE, Population Method)

```python
import pyexocross as px

px.stick_spectra(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exoatom_stick_nlte.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['configuration', 'Multiple', 'parity'],
    qnsformat_list=['%20s', '%10s', '%2s'],
    nlte_method='P',
    nlte_path='/path/to/Databases/ExoAtom/Ar/NIST/Ar_Ids.csv',
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=115400,
    abs_emi='Ab',
    qns_filter={
        'configuration': [],
        'Multiple': [],
        'parity': [],
    },
    plot=True,
    plot_method='log',
    plot_wn_wl='WN', 
    plot_unit='cm-1',
    limit_yaxis=1e-30,
)
```

### Cross Sections (Non-LTE, Population Method)

```python
import pyexocross as px

px.cross_sections(
    database='ExoAtom',
    atom='Ar',
    dataset='NIST',
    species_id=601,
    read_path='/path/to/Databases/ExoAtom/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/exoatom_xsec_nlte.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['configuration', 'Multiple', 'parity'],
    qnsformat_list=['%20s', '%10s', '%2s'],
    nlte_method='P',
    nlte_path='/path/to/Databases/ExoAtom/Ar/NIST/Ar_Ids.csv',
    pressures=[1.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=115400,
    bin_size=0.1,
    profile='SciPyVoigt',
    broadeners=['Default'],
    ratios=[1.0],
    abs_emi='Ab',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN', 
    plot_unit='cm-1',
    limit_yaxis=1e-30,
)
```

---

## HITRAN Examples

### Conversion (HITRAN -> ExoMol)

```python
import pyexocross as px

px.conversion(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran.log',
    chunk_size=100000,
    qnslabel_list=['J', 'X', 'Omega', 'v1', 'Sym', 'F'],
    qnsformat_list=['%5.1f', '%2s', '%3s', '%2d', '%1s', '%5s'],
    conversion_format=2,
    conversion_min_freq=0,
    conversion_max_freq=63000,
    conversion_unc=None,           
    conversion_threshold=None,     
    global_qn_label_list=['X', 'Omega', 'v1'],
    global_qn_format_list=['%2s', '%3s', '%2d'],
    local_qn_label_list=['Br', 'Sym', 'F'],
    local_qn_format_list=['%2s', '%1s', '%5s'],
)
```

### Partition Functions

```python
import pyexocross as px

px.partition_functions(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran_pf.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['J', 'X', 'Omega', 'v1', 'Sym', 'F'],
    qnsformat_list=['%5.1f', '%2s', '%3s', '%2d', '%1s', '%5s'],
    ntemp=1,
    tmax=5000,
)
```

### Specific Heats

```python
import pyexocross as px

px.specific_heats(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran_cp.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['J', 'X', 'Omega', 'v1', 'Sym', 'F'],
    qnsformat_list=['%5.1f', '%2s', '%3s', '%2d', '%1s', '%5s'],
    ntemp=1,
    tmax=5000,
)
```

### Cooling Functions

```python
import pyexocross as px

px.cooling_functions(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran_cf.log',
    chunk_size=100000,
    ntemp=1,
    tmax=5000,
)
```

### Lifetimes

```python
import pyexocross as px

px.lifetimes(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran_lifetime.log',
    ncputrans=4,    
    ncpufiles=1,
    chunk_size=100000,
    qnslabel_list=['J', 'X', 'Omega', 'v1', 'Sym', 'F'],
    qnsformat_list=['%5.1f', '%2s', '%3s', '%2d', '%1s', '%5s'],
    compress=False,
)
```

### Oscillator Strengths

```python
import pyexocross as px

px.oscillator_strengths(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran_os.log',
    chunk_size=100000,
    gf_or_f='f',
    plot=True,
    plot_method='log',
    plot_wn_wl='WN', 
    plot_unit='cm-1',
    limit_yaxis=1e-30,
)
```

### Stick Spectra (Non-LTE, Treanor Method)

```python
import pyexocross as px

px.stick_spectra(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran_stick_nlte.log',
    chunk_size=100000,
    qnslabel_list=['J', 'X', 'Omega', 'v1', 'Sym', 'F'],
    qnsformat_list=['%5.1f', '%2s', '%3s', '%2d', '%1s', '%5s'],
    nlte_method='T',
    tvib_list=[1000, 2000],
    trot_list=[300, 400],
    vib_label=['v1', 'X'],
    rot_label=['J'],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    abs_emi='Ab',
    unc_filter=0.01,
    threshold=1e-35,
    plot=True,
    plot_method='log',
    plot_wn_wl='WN', 
    plot_unit='cm-1',
    limit_yaxis=1e-30,
)
```

### Cross Sections (with Air + Self Broadening)

```python
import pyexocross as px

px.cross_sections(
    database='HITRAN',
    molecule='NO',
    isotopologue='14N-16O',
    dataset='NO-HITRAN',
    species_id=81,
    read_path='/path/to/Databases/HITRAN/',
    save_path='/path/to/output/',
    logs_path='/path/to/output/log/hitran_xsec.log',
    chunk_size=100000,
    qnslabel_list=['J', 'X', 'Omega', 'v1', 'Sym', 'F'],
    qnsformat_list=['%5.1f', '%2s', '%3s', '%2d', '%1s', '%5s'],
    temperatures=[1000, 2000],
    pressures=[0.5, 1.0],
    wn_wl='WN',
    wn_wl_unit='cm-1',
    min_range=0,
    max_range=30000,
    bin_size=0.1,
    profile='SciPyVoigt',
    broadeners=['Air', 'Self'],
    ratios=[0.7, 0.3],
    cutoff=25.0,
    abs_emi='Ab',
    plot=True,
    plot_method='log',
    plot_unit='cm-1',
    limit_yaxis=1e-30,
)
```

---

## Batch Processing Example

Run multiple species or temperatures in a loop:

```python
import pyexocross as px

species_list = [
    {'molecule': 'NO', 'isotopologue': '14N-16O', 'dataset': 'XABC', 'species_id': 81},
    {'molecule': 'NO', 'isotopologue': '15N-16O', 'dataset': 'XABC', 'species_id': 82},
]

COMMON = dict(
    database='ExoMol',
    read_path='/home/jingxin/LHD/Program/Databases/ExoMol/',
    save_path='/home/jingxin/LHD/Program/Data/pyexocross/',
    logs_path='/home/jingxin/LHD/Program/Data/pyexocross/log/test_api_exomol.log',
)

COMPUTE_PARAMS = dict(
    ncputrans=4,                
    ncpufiles=1,               
    chunk_size=100000,         
)

QN_PARAMS = dict(
    qnslabel_list=['par', 'e/f', 'eS', 'v', 'Lambda', 'Sigma', 'Omega'],
    qnsformat_list=['%1s', '%1s', '%13s', '%3d', '%2d', '%7.1f', '%7.1f'],
)

RANGE_PARAMS = dict( 
    wn_wl='WN',                 
    wn_wl_unit='cm-1',          
    min_range=0,               
    max_range=30000,            
    abs_emi='Ab',              
    unc_filter=0.01,          
    threshold=1e-30,           
)

LINE_PROFILE = dict(
    profile='SciPyVoigt',
    broadeners=['Default'], 
    ratios=[1.0],           
    alpha_hwhm=None,         
    gamma_hwhm=None,   
)

PLOT = dict(
    plot=True,             
    plot_method='log',     
    plot_wn_wl='WN',       
    plot_unit='cm-1',       
    limit_yaxis=1e-30,      
)

for species in species_list:
    px.cross_sections(
        **COMMON,
        **COMPUTE_PARAMS,
        **QN_PARAMS,
        **RANGE_PARAMS,
        **LINE_PROFILE,
        **PLOT,
        temperatures=[1000, 3000, 5000],
        pressures=[1.0, 10.0],
        bin_size=0.1,           
        predissociation=False,  
        cutoff=25.0,            
    )
```
