# Basic information

Basic information includes data source, file path and functions.

## Data source

In data source section, please provide:

1. The name of database, molecule, isotopologue and dataset.
2. The molecule ID and isotopologue ID.

### Molecule database

**For ExoMol database**

The name of `Database`, `Molecule`, `Isotopologue`, and `Dataset` are necessary.

The molecule and isotopologue ID `SpeciesID` can be set as `0` or any other integers.

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                H2O
Isotopologue                            1H2-16O
Dataset                                 POKAZATEL
SpeciesID                               11
```

**For ExoMol database**

The name of `Database`, `Molecule`, and `Isotopologue` are necessary.

The molecule and isotopologue ID `SpeciesID` can be set as `0` or any other integers.

*Example*

```bash
# Data source #
Database                                ExoMolHR
Molecule                                NO
Isotopologue                            14N-16O
SpeciesID                               81
```

**For HITRAN and HITEMP databases**

The `Database` name, `Molecule` name, and the molecule and isotopologue ID `SpeciesID` are necessary.

The first two digits of `SpeciesID` are molcule ID and the third digit is isotopologue ID and there is no blank between molecule ID and isotopologue ID. `SpeciesID` can be found from [HITRANOnline Isotopologue Metadata](https://hitran.org/docs/iso-meta/).

The name of `Isotopologue` and `Dataset` can be set as 'none' or any other strings.

*Example*

```bash
# Data source #
Database                                HITRAN
Molecule                                CO2
Isotopologue                            none
Dataset                                 none
SpeciesID                               21
```

```bash
# Data source #
Database                                HITRAN
Molecule                                C2H2
Isotopologue                            abcd
Dataset                                 EfGh
SpeciesID                               261
```

### Atom and ion database

**For ExoAtom database**

The name of `Database`, `Atom`, and `Dataset` are necessary.

Only provide two datasets: `NIST` and `Kurucz`.

*Example*

```bash
# Data source #
Database                                ExoAtom
Atom                                    Ar
Dataset                                 NIST
```

```bash
# Data source #
Database                                ExoAtom
Atom                                    Al
Dataset                                 Kurucz
```

## File path

File path section records the file path for both reading and saving.

`ReadPath` and `SavePath` are folder paths and should start with `/`, but don't need to end with `/`.

`LogFilePath` is the file path of the log file, not the folder path.

&#x2705; /aaa/bbb/ccc

&#x2705; /aaa/bbb/ccc/

&#x2705; /aaa/bbb/ccc/ddd.log

**For ExoMol database**

`ReadPath` is the folder path of input files when the input line list files path is stored in the following format.

<font color=Orange>`ReadPath`</font>/<font color=Pink>`Molecule`</font>/<font color=SkyBlue>`Isotopologue`</font>/<font color=YellowGreen>`Dataset`</font>/<font color=SkyBlue>`Isotopologue`</font><font color=Brown>__</font><font color=YellowGreen>`Dataset`</font><font color=Brown>.states.bz2</font>

<font color=Orange>/mnt/data/exomol/exomol3_data</font>/<font color=Pink>MgH</font>/<font color=SkyBlue>24Mg-1H</font>/<font color=YellowGreen>XAB</font>/<font color=SkyBlue>24Mg-1H</font><font color=Brown>__</font><font color=YellowGreen>XAB</font><font color=Brown>.states.bz2</font></font>

```
└── exomol3_data
           ├── C2H2
           ├── CO2
           │     ├── 12C-16O2
           │     ├── 13C-16O2
           │     ├── ...
           │     ├── 12C-16O2__air.broad
           │     └── 12C-16O2__self.broad
           ├── MgH
           │     ├── 24Mg-1H
           │     │       ├── Yadin
           │     │       └── XAB
           │     │            ├── 24Mg-1H__XAB.def
           │     │            ├── 24Mg-1H__XAB.def.json
           │     │            ├── 24Mg-1H__XAB.pf
           │     │            ├── 24Mg-1H__XAB.states.bz2
           │     │            └── 24Mg-1H__XAB.trans.bz2
           │     ├── 25Mg-1H
           │     │       ├── Yadin
           │     │       └── XAB
           │     │            ├── 25Mg-1H__XAB.def
           │     │            ├── 25Mg-1H__XAB.def.json
           │     │            ├── 25Mg-1H__XAB.pf
           │     │            ├── 25Mg-1H__XAB.states.bz2
           │     │            └── 25Mg-1H__XAB.trans.bz2
           │     └── 26Mg-1H
           ├── ...
           │
```

`SavePath` is the folder path for saving all results obtained by the PyExoCross program.

`LogFilePath` is the file path of the log file, the program can record the log output automatically.

*Example*

```bash
# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/MgH_ExoMol.log
```

**For ExoMolHR database**

`ReadPath` is the folder path of input files when the input line list files path is stored in the following format.

For ExoMolHR line list files, the standard filename formats are: 

<font color=Orange>`ReadPath`</font>/<font color=Pink>`Molecule`</font>/<font color=SkyBlue>`Isotopologue`</font>/<font color=Brown></font><font color=SkyBlue>date__`Isotopologue`__T</font><font color=Brown>.csv</font>

<font color=Orange>/mnt/data/exomolhr/exomolhr_results</font>/<font color=Pink>NO</font>/<font color=SkyBlue>14N-16O</font>/<font color=Brown>20260311080614__</font><font color=SkyBlue>14N-16O</font><font color=Brown>__1000K.csv</font></font>

<font color=Orange>`ReadPath`</font>/<font color=Pink>`Molecule`</font>/<font color=SkyBlue>`Isotopologue`</font>/<font color=SkyBlue>`Isotopologue`</font><font color=Brown>__</font><font color=YellowGreen>`Dataset`</font><font color=Brown>.pf</font>

<font color=Orange>/mnt/data/exomol/exomol3_data</font>/<font color=Pink>NO</font>/<font color=SkyBlue>14N-16O</font>/<font color=SkyBlue>14N-16O</font><font color=Brown>__</font><font color=YellowGreen>XABC</font><font color=Brown>.pf</font></font>

However, users can also rename the CSV filenames in any formats as long as `Isotopologue` is included in the filename.

&#x2705; date__`Isotopologue`__T.csv

&#x2705; `Molecule`__`Isotopologue`.csv

&#x2705; `Molecule`__`Isotopologue`__T.csv

```
└── exomolhr_results
           ├── C2H2
           ├── CO2
           │     ├── 12C-16O2
           │     ├── 13C-16O2
           │     ├── ...
           │     ├── 12C-16O2__air.broad
           │     └── 12C-16O2__self.broad
           ├── NO
           │     └── 14N-16O
           │             ├── 14N-16O__XABC.pf
           │             └── 20260311080614__14N-16O__1000K.csv
           ├── MgH
           │     ├── 24Mg-1H
           │     │       ├── 24Mg-1H__XAB.pf
           │     │       └── 20260311080614__24Mg-1H__1000K.csv
           │     ├── 25Mg-1H
           │     │       ├── 25Mg-1H__XAB.pf
           │     │       └── MgH__25Mg-1H__1000K.csv
           │     └── 26Mg-1H    
           │             ├── 26Mg-1H__XAB.pf
           │             └── MgH__26Mg-1H.csv     
           ├── ...
           │
```

`SavePath` is the folder path for saving all results obtained by the PyExoCross program.

`LogFilePath` is the file path of the log file, the program can record the log output automatically.

*Example*

```bash
# File path #
ReadPath                                /mnt/data/exomolhr/exomolhr_results/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_ExoMolHR.log
```

**For ExoAtom database**

`ReadPath` is the folder path of input files when the input line list files path is stored in the following format.

<font color=Orange>`ReadPath`</font>/<font color=SkyBlue>`Atom`</font>/<font color=YellowGreen>`Dataset`</font>/<font color=SkyBlue>`Atom`</font><font color=Brown>__</font><font color=YellowGreen>`Dataset`</font><font color=Brown>.states</font>

<font color=Orange>/mnt/data/exoatom/exoatom_data</font>/<font color=SkyBlue>Li</font>/<font color=YellowGreen>NIST</font>/<font color=SkyBlue>Li</font><font color=Brown>__</font><font color=YellowGreen>NIST</font><font color=Brown>.states</font></font>

`SavePath` is the folder path for saving all results obtained by the PyExoCross program.

`LogFilePath` is the file path of the log file, the program can record the log output automatically.

```
└── exoatom_data
           ├── Al
           ├── Al_p
           ├── Li
           │    ├── NIST
           │    │     ├── Li_p__NIST.adef.json
           │    │     ├── Li_p__NIST.pf
           │    │     ├── Li_p__NIST.states
           │    │     └── Li_p__NIST.trans
           │    └── Kurucz
           │          ├── Li_p__Kurucz.adef.json
           │          ├── Li_p__Kurucz.pf
           │          ├── Li_p__Kurucz.states
           │          └── Li_p__Kurucz.trans
           ├── Li_p
           │    ├── NIST
           │    │     ├── Li__NIST.adef.json
           │    │     ├── Li__NIST.pf
           │    │     ├── Li__NIST.states
           │    │     └── Li__NIST.trans
           │    └── Kurucz
           │          ├── Li__Kurucz.adef.json
           │          ├── Li__Kurucz.pf
           │          ├── Li__Kurucz.states
           │          └── Li__Kurucz.trans   
           ├── ...
           │
```

*Example*

```bash
# File path #
ReadPath                                /mnt/data/exoatom/exoatom_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/Li_NIST.log
```

**For HITRAN and HITEMP databases**

`ReadPath` is the folder path of input files when the input line list files path is stored in the following format.

<font color=Orange>`ReadPath`</font>/<font color=SkyBlue>`Molecule`</font>/<font color=YellowGreen>`Isotopologue`</font>/<font color=SkyBlue>`Molecule`</font><font color=Brown>__</font><font color=YellowGreen>`Isotopologue`</font><font color=Brown>.par</font>

<font color=Orange>/path/HITRAN</font>/<font color=SkyBlue>NO</font>/<font color=YellowGreen>14N-16O</font>/<font color=SkyBlue>NO</font><font color=Brown>__</font><font color=YellowGreen>14N-16O</font><font color=Brown>.par</font></font>

`SavePath` is the folder path for saving all results obtained by the PyExoCross program.

`LogFilePath` is the file path of the log file, the program can record the log output automatically.

```
└── HITRAN
           ├── H2O
           ├── CO2
           ├── NO
           │    ├── 14N-16O
           │    │     ├── NO__14N-16O.par
           │    │     └── NO__14N-16O.pf / NO__14N-16O.txt
           │    ├── 15N-16O
           │    │     ├── NO__15N-16O.par
           │    │     └── NO__15N-16O.pf / NO__15N-16O.txt
           │    └── 14N-18O
           │          ├── NO__14N-18O.par
           │          └── NO__14N-18O.pf / NO__14N-18O.txt
           ├── NO_p
           │    └── 14N-16O_p
           │          ├── NO_p__14N-16O_p.par
           │          └── NO_p__14N-16O_p.pf / NO_p__14N-16O_p.txt
           ├── ...
           │
```

*Example*

```bash
# File path #
ReadPath                                /home/jingxin/data/HITRAN/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/CO2_HITRAN.log
```

## Functions

In current version, *PyExoCross* can convert data format between the ExoMol (ExoMol, ExoMolHR, and ExoAtom databases) and the HITRAN (HITRAN and HITEMP database) formats.

*PyExoCross* also implements the computations of other useful functions including partition functions, specific heats, cooling functions, radiative lifetimes, oscillator strengths, LTE and non-LTE stick spectra and cross sections.

*PyExoCross* can only compute stick spectra and cross sections for ExoMolHR database.

If users use `.inp` input files:

*PyExoCross* can calculate partition functions, cooling functions, oscillator strengths, LTE and non-LTE stick spectra and cross sections directly for data from the HITRAN format databases. Specific heats and lifetimes still require HITRAN data to be converted to the ExoMol format first and then treated as ExoMol data in *PyExoCross*. 

If users use Python package:

*PyExoCross* can convert data format automatically if required.


Use this function or not:

`0` means no;

`1` means yes.

If the value of a function's second column is `0`, then there is no need to do any changes in this function's own section, the program won't process data with this function.

*Example*

```bash
# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           1
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           1
```

## Cores and chunks

Please provide the number of cores `NCPU` and the size of chunks `ChunkSize` of the quantum numbers when you use *PyExoCross* uses multiprocessing.

The program will run on different cores together. 

`NCPUtrans`: The number of cores for processing each transitions file. 

`NCPUfiles`: The number of transitions files for processing at the same time. 

`ChunkSize`: The program splits each transitions file to many chunks when reading and calculating. `ChunkSize` is the size of each chunk.

`RunMode`: Choose to run the program in CPU or GPU mode.

`GPUBackend`: GPU backend selection (only used when `RunMode=GPU`):

- `'AUTO'` (recommended): `PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback`
- `'CUDA'`: `PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback`
- `'PyTorch-CUDA'`: NVIDIA PyTorch CUDA only; otherwise CPU fallback
- `'CuPy-CUDA'`: NVIDIA CuPy CUDA only; otherwise CPU fallback
- `'MPS'`: Apple Metal (MPS) only; otherwise CPU fallback

**Note**

- MPS uses float32 kernels, so results can differ slightly from CPU/CUDA float64.
- If no compatible GPU backend is available, PyExoCross falls back to CPU automatically.
- GPU acceleration is available for cooling functions, stick spectra, cross sections, and stick spectra + cross sections.


`GPUBatchLines`: The number of lines to process in each batch on the GPU. Default value is 8192.

`GPUBatchGrid`: The number of grids to process in each batch on the GPU. Default value is 256.

**Note**

Some suggestions on setting the number of `NCPUtrans` and `NCPUfiles`.

1. `NCPUtrans` $*$ `NCPUfiles` $\leq$ Your cores number
2. Some cases (depend on `.trans` files):
   
   | `.trans` files size | `.trans` files number | `NCPUtrans` ? `NCPUfiles`         |
   | :-----------------: | :-------------------: | :-------------------------------: |
   | Very small          |  Not small            | `NCPUtrans` $=1<$ `NCPUfiles`     |
   | Not large           |  Very large           | `NCPUtrans` $<$ `NCPUfiles`       |
   | Not large           |  Not large            | `NCPUtrans` $\approx$ `NCPUfiles` |
   | Very large          |  Small                | `NCPUtrans` $>1=$ `NCPUfiles`     |

*Example*

```bash
# Cores and chunks #
NCPUtrans                               2
NCPUfiles                               4
ChunkSize                               500000
RunMode                                 CPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default): PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback
GPUBatchLines                           8192                      # GPU line-batch size (only used when RunMode=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when RunMode=GPU)
```


```bash
# Cores and chunks #
NCPUtrans                               8
NCPUfiles                               1
ChunkSize                               1000000
RunMode                                 GPU                       # CPU(default) or GPU
GPUBackend                              AUTO                      # AUTO(default): PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback
GPUBatchLines                           8192                      # GPU line-batch size (only used when RunMode=GPU)
GPUBatchGrid                            256                       # GPU grid-batch size (only used when RunMode=GPU)
```

```bash
# CUDA policy (NVIDIA)
# Priority: PyTorch-CUDA -> CuPy-CUDA -> MPS -> CPU fallback
RunMode                                 GPU
GPUBackend                              CUDA
GPUBatchLines                           8192
GPUBatchGrid                            256
```

```bash
# Force PyTorch CUDA only
RunMode                                 GPU
GPUBackend                              PyTorch-CUDA
GPUBatchLines                           8192
GPUBatchGrid                            256
```

```bash
# Force CuPy CUDA only
RunMode                                 GPU
GPUBackend                              CuPy-CUDA
GPUBatchLines                           8192
GPUBatchGrid                            256
```

```bash
# Force MPS (Apple Silicon)
RunMode                                 GPU
GPUBackend                              MPS
GPUBatchLines                           8192
GPUBatchGrid                            256
```

## Quantum numbers

Please provide the labels `QNslabel` and formats `QNsformat` of the quantum numbers when you use *PyExoCross* to convert data format, calculate stick spectra or cross sections if you need the quantum filter.

* The definition files `.def`, `.def.json`, and `.adef.json` of ExoMol and ExoAtom databases (available at [exomol.com](https://www.exomol.com/)) provides the labels and formats of the quantum numbers for each species for reference.
* HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.

**Note**

You can define the quantum number column name by yourself, but please make sure it has letters without any blanks. \
e.g. 'c1', 'c2', 'v1', 'v2', 'electronicState', 'electronic_state', '1v', '2v', 'M/E/C'. \
Wrong format of the quantum number column nams: '1', '2', 'electronic state'.

*Example*

```bash
# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %1d      %7.1f    %7.1f
```
