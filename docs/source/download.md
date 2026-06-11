# Database Downloading

PyExoCross provides a built-in Python API for downloading molecular and atomic data from external databases, including **ExoMol**, **ExoAtom**, **ExoMolHR**, and **HITRAN**.

---

## Python API Usage

All database downloads are handled through the main function `pyexocross.download()`.

```python
import pyexocross as px

px.download(database, file_path=None, **kwargs)
```

### Common Parameters

| Parameter | Type | Required | Description |
|---|---|---|---|
| `database` | `str` | Yes | Database name: `'ExoMol'`, `'ExoAtom'`, `'ExoMolHR'`, or `'HITRAN'` |
| `file_path` | `str` | Yes | Folder where the downloaded files will be saved, can also use `save_path` |
| `species_info` | `dict` | Yes | Molecule/isotopologue or atom configuration dict |
| `download` | `bool` | yes | If `True` (default), download the files. If `False`, only generate the URL list file. |
| `wn_range` | `list` | Optional | Filter transitions by wavenumber range (cm⁻¹). |
| `T` | `int` | Optional | Temperature in Kelvin. |
| `threshold` | `float` | Optional | Minimum line intensity threshold in cm/molecule. |
| `dataset` | `str` | Optional | The dataset to use. Can be either `'NIST'` or `'Kurucz'`. |

### Species Configuration Details (`species_info`)

The dictionary structure for `species_info` allows you to customize the files to be downloaded. You can specify parameters such as:

*   **`wn_range`** (list of two floats: `[wn_min, wn_max]`): Filter transition files by wavenumber range (in cm⁻¹). 
    *   *ExoMol*: Used to download only the segmented transition files (`.trans.bz2`) that fall completely within the requested range.
    *   *HITRAN* / *ExoMolHR*: Used to query and download the line transitions matching this specific range.
*   **`T`** or **`temperature`** (int): Temperature (in Kelvin) for fetching temperature-dependent high-resolution datasets.
    *   *ExoMolHR*: Required for online query-based CSV download.
*   **`threshold`** or **`Smin`** (float): The minimum line intensity threshold (in cm/molecule).
    *   *ExoMolHR*: Required for online query-based CSV download to filter weak lines.
*   **`dataset`** (str): The dataset to use. Required for ExoAtom downloads. 
    *   *ExoAtom*: Can be either `'NIST'` (critically evaluated and recommended for accuracy) or `'Kurucz'` (recommended for completeness).

| Database | Required Parameters | Optional Parameters |
| :--- | :--- | :--- |
| `'ExoMol'` | `database`, `file_path`, `species_info`, `download` | `wn_range` |
| `'ExoMolHR'` | `database`, `file_path`, `species_info`, `download` | `T`, `threshold`, `wn_range` |
| `'ExoAtom'`| `database`, `file_path`, `species_info`, `download` | `dataset` |
| `'HITRAN'` | `database`, `file_path`, `species_info`, `download` | `wn_range` |

---

## 1. ExoMol Database (`database='ExoMol'`)

Downloads line lists, states, transition files, and partition functions from the [ExoMol database](www.exomol.com).

### Example

```python
px.download(
    database='ExoMol',
    file_path='/path/to/ExoMol/',
    species_info={
        'MgH': {
            '24Mg-1H': {'wn_range': None},
            '25Mg-1H': {'wn_range': None},
        },
        'H2O': {
            '1H2-16O': {'wn_range': [41000, 41200]},
        },
    },
    download=True,
)
```

---

## 2. ExoMolHR Database (`database='ExoMolHR'`)

Downloads high-resolution line lists and partition functions. It queries the server for a dynamically calculated CSV file based on your specified temperature (`T`), wavenumber range (`wn_range`), and intensity threshold (`threshold`).

### Example

```python
px.download(
    database='ExoMolHR',
    file_path='/path/to/ExoMolHR/',
    species_info={
        'MgH': {
            '24Mg-1H': {'T': 1000, 'wn_range': [0, 500], 'threshold': 1e-30},
            '25Mg-1H': None,
        },
        'AlH': {
            '27Al-1H': {'T': 500, 'wn_range': [0, 500], 'threshold': 1e-30}
        },
    },
    download=True,
)
```

---

## 3. ExoAtom Database (`database='ExoAtom'`)

Downloads atomic transitions and energy levels from the NIST or Kurucz databases using the ExoAtom configuration.

### Example

```python
px.download(
    database='ExoAtom',
    file_path='/path/to/ExoAtom/',
    species_info={
        'He': {'dataset': 'NIST'},
        'He_p': {'3He_p': {'dataset': 'NIST'}},
        'Ar_p': {'dataset': 'Kurucz'},
    },
    download=True,
)
```

---

## 4. HITRAN Database (`database='HITRAN'`)

Downloads line-by-line spectroscopic data and partition functions directly from the [HITRAN database](https://hitran.org) via the HITRANonline API.

### Example

```python
px.download(
    database='HITRAN',
    file_path='/path/to/HITRAN/',
    species_info={
        'NO': {
            '14N-16O': {'wn_range': [0, 100]},
            '15N-16O': {'wn_range': [100, 150]},
        },
        'H2O': {
            '1H2-16O': {'wn_range': [100, 110]},
        },
    },
    download=True,
)
```

---

## Generating URLs Only

If you only want to generate the list of URLs to be downloaded later, set `download=False`:

```python
px.download(
    database='ExoMol',
    file_path='/path/to/ExoMol/',
    species_info={'MgH': {'24Mg-1H': None}},
    download=False,
)
```

This will save the URLs to `/path/to/ExoMol/url/exomol__urls.txt` and print them to the console.
