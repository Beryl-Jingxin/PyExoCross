# Interactive results

This guide explains how to keep PyExoCross results in Python for use in
Jupyter, Spyder, or a Python script. Start with the short example below, then
use the section that matches your calculation.

## 1. Choose an output mode

All calculation functions support the same three modes:

| Mode | Save result files | Return a `Result` | Save plots |
|---|---:|---:|---:|
| `output='files'` | Yes | No | When `plot=True` |
| `output='memory'` | No | Yes | No |
| `output='both'` | Yes | Yes | When `plot=True` |

`files` is the default and preserves the original PyExoCross behaviour.
Use `memory` while changing parameters and comparing results interactively.
Use `both` when the same run should also create permanent files.

```python
result = px.cross_sections(
    data=data,
    temperatures=[300],
    pressures=[1.0],
    output='memory',
    log='none',
)
```

***Important***

`px.conversion()` is a format transformation. It supports `output='files'`
only and returns `None`.


## 2. Load once, calculate many times

`px.load()` returns a reusable `LoadedData` object. It does not return a
calculation `Result`.

```python
import pyexocross as px

data = px.load(
    database='ExoMol',
    molecule='MgH',
    isotopologue='24Mg-1H',
    dataset='XAB',
    read_path='/path/to/ExoMol/',
    min_range=0,
    max_range=10000,
    cache='auto',
    log='none',
)

result = px.cross_sections(
    data=data,
    temperatures=[300, 500],
    pressures=[0.1, 1.0],
    min_range=0,
    max_range=100,
    bin_size=0.1,
    profile='Gaussian',
    output='memory',
    log='none',
)
```

The calculation range must remain inside the range loaded by `px.load()`.
Lifetimes, cooling functions, and oscillator strengths automatically load all
transitions when required.

You can also call any calculation directly without `px.load()`. In that case,
provide `database`, species information, and `read_path` to the calculation.

## 3. What is inside `Result`?

Type `result.` and press Tab in IPython or Spyder, or run:

```python
result.help()
result.summary()
```

The public fields are:

| Field | Meaning |
|---|---|
| `result.data` | Calculated arrays or tables grouped by physical condition |
| `result.coords` | Temperature, state, wavenumber, or wavelength coordinates |
| `result.units` | Coordinate and data units |
| `result.conditions` | Available `T`, `P`, `Tvib`, and `Trot` combinations |
| `result.params` | Resolved parameters relevant to this calculation |
| `result.timing` | Load, prepare, calculate, save, and total times |
| `result.select()` | Select one calculated array or table |

Normally, use `select()` instead of reading `result.data` directly.

## 4. Select results by calculation type

### 4.1 One-table calculations

Partition functions, specific heats, cooling functions, lifetimes, and
oscillator strengths produce one table. Select it without conditions:

```python
result = px.partition_functions(
    data=data,
    ntemp=100,
    tmax=5000,
    output='memory',
    log='none',
)

table = result.select()
print(table.head())
print(result.coords['temperature'][:5])
```

The same access pattern applies to:

```python
specific = px.specific_heats(
    data=data, ntemp=100, tmax=5000, output='memory', log='none'
).select()

cooling = px.cooling_functions(
    data=data, ntemp=100, tmax=5000, output='memory', log='none'
).select()

lifetime = px.lifetimes(
    data=data, output='memory', log='none'
).select()

oscillator = px.oscillator_strengths(
    data=data, gf_or_f='f', output='memory', log='none'
).select()
```

### 4.2 Stick spectra

Stick spectra depend on temperature but not pressure:

```python
result = px.stick_spectra(
    data=data,
    temperatures=[300, 1000],
    output='memory',
    log='none',
)

stick_300 = result.select(T=300)
stick_1000 = result.select(T=1000)
```

### 4.3 Cross sections

Cross sections are selected by temperature and pressure:

```python
result = px.cross_sections(
    data=data,
    temperatures=[300, 1000],
    pressures=[0.1, 1.0],
    output='memory',
    log='none',
)

xsec = result.select(T=300, P=1.0)
```

If a selection matches more than one result, PyExoCross asks for more
conditions. Use `squeeze=False` when all matches are wanted:

```python
all_pressures_at_300 = result.select(T=300, squeeze=False)
```

This returns a mapping from `Condition` to calculated data.

### 4.4 Combined stick spectra and cross sections

A combined calculation contains two products. Always specify `product`:

```python
result = px.stick_spectra_cross_section(
    data=data,
    temperatures=[300],
    pressures=[1.0],
    output='memory',
    log='none',
)

stick = result.select(T=300, product='stick_spectra')
xsec = result.select(T=300, P=1.0, product='cross_section')
```

Stick spectra do not use pressure. Their condition therefore has `P=None`.

Combined units are also separated by product:

```python
result.units['stick_spectra']
# {'wavenumber': 'cm-1', 'data': 'cm/molecule'}

result.units['cross_section']
# {'wavenumber': 'cm-1', 'data': 'cm2/molecule'}
```

### 4.5 Non-LTE results

For the two-temperature method, select with `Tvib` and `Trot`:

```python
result = px.stick_spectra(
    data=data,
    nlte_method='T',
    tvib_list=[3000],
    trot_list=[296],
    vib_label=['v', 'ElecState'],
    rot_label=['J', 'e/f'],
    output='memory',
    log='none',
)

stick = result.select(Tvib=3000, Trot=296)
```

For combined Non-LTE output, add `product='stick_spectra'` or
`product='cross_section'`.

## 5. Coordinates and units

For a cross section with one shared grid:

```python
grid = result.coords['wavenumber']
unit = result.units['wavenumber']
data_unit = result.units['data']
```

When conditions use different grids, the coordinate value is a mapping from
`Condition` to arrays. Inspect it with:

```python
print(result.coords)
print(result.conditions)
```

Use the condition matching the data selected from `result.select()`.

## 6. Calculation parameters

`result.params` contains only parameters relevant to the selected function and
active calculation mode. For example, a partition-function result does not
show pressure, profile, or broadening settings. LTE results do not show unused
Non-LTE parameters.

```python
print(result.params)                 # readable complete list
profile = result.params['profile']   # one parameter
params = result.params.to_dict()     # ordinary dictionary
```

Values are resolved values used internally. For example, `Gaussian` may be
shown as `GAUSSIAN`.

## 7. Timing

```python
print(result.timing)
calculate_wall = result.timing['calculate']['system']
save_wall = result.timing['save']['system']
```

| Timing key | Included work |
|---|---|
| `load` | Reading source data for this calculation call |
| `prepare` | Filtering and preparing arrays/tables |
| `calculate` | Physical calculation; excludes result-file saving and plotting |
| `save` | Writing result files and plots |
| `total` | Complete calculation call |

When preloaded `data` is reused, `load` may be zero because loading happened in
the earlier `px.load()` call.

## 8. Help and discovery

```python
px.help()                         # overview
px.help('cross_sections')         # one function
px.help(px.cross_sections)        # equivalent function-object form
px.help('result')                 # Result fields and selection
result.help()                     # examples for this result
help(px.cross_sections)           # Python built-in help
```

In IPython and Jupyter, these are also available:

```python
px.cross_sections?
px.cross_sections??
```

## 9. Logging and terminal output

Logging and terminal output are independent of `output`:

| Option | Behaviour |
|---|---|
| `log='auto'` | Write a log when `logs_path` or `.inp` `LogFilePath` is valid |
| `log='file'` | Require a valid log path |
| `log='none'` | Do not create a log file |
| `verbose=True` | Show normal terminal output and progress |
| `verbose=False` | Hide normal terminal output; an enabled log still receives it |

A quiet interactive calculation is therefore:

```python
result = px.cross_sections(
    data=data,
    temperatures=[300],
    pressures=[1.0],
    output='memory',
    log='none',
    verbose=False,
)
```

## 10. Common mistakes

| Problem | Correct use |
|---|---|
| `result = px.cross_sections` stores the function | Call it: `result = px.cross_sections(...)` |
| `result` is `None` | Use `output='memory'` or `output='both'` |
| More than one condition matches | Add `T`, `P`, `Tvib`, or `Trot`; or use `squeeze=False` |
| Combined selection says product is required | Add `product='stick_spectra'` or `product='cross_section'` |
| Stick selection finds no pressure match | Remove `P`; stick spectra do not use pressure |
| Conversion rejects memory output | Use `output='files'` |

## 11. Complete notebooks

The repository contains two executable tutorials based on the ExoMol MgH
parameters in `test/test_api_exomol.py`:

- `test/tutorial_result_api.ipynb`: load, conversion, every calculation
  function, Result inspection, selection, parameters, units, and timing.
- `test/tutorial_help_api.ipynb`: package help, function help, built-in Python
  help, Result help, and common mistakes.
