# Basic information

Basic information includes data source, file path and functions.

## Data source

In data source section, please provide:

1. The name of database, molecule, isotopologue and dataset.
2. The molecule ID and isotopologue ID.

**For ExoMol**

The name of `Database`, `Molecule`, `Isotopologue` and `Dataset` are necessary.

The molecule and isotopologue ID `MolIsoID` can be set as `0` or any other integers.

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
MolIsoID                                501
```

**For HITRAN**

The `Database` name, `Molecule` name and the molecule and isotopologue ID `MolIsoID` are necessary. The first two digits of `MolIsoID` are molcule ID and the third digit is isotopologue ID and there is no blank between molecule ID and isotopologue ID.

The name of `Isotopologue` and `Dataset` can be set as 'none' or any other strings.

*Example*

```bash
# Data source #
Database                                HITRAN
Molecule                                CO2
Isotopologue                            none
Dataset                                 none
MolIsoID                                21
```

```bash
# Data source #
Database                                HITRAN
Molecule                                C2H2
Isotopologue                            abcd
Dataset                                 EfGh
MolIsoID                                261
```

## File path

File path section records the file path for both reading and saving.

**For ExoMol**

`ReadPath` is the input files' folder path when the input line list files path is stored in the following format.

```
/FolderPath/molecule/iso-slug/dataset/iso-slug__dataset.states.bz2

/mnt/data/exomol/exomol3_data/MgH/24Mg-1H/XAB/24Mg-1H__XAB.states.bz2
```

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
           │     │            ├── 24Mg-1H__XAB.pf
           │     │            ├── 24Mg-1H__XAB.states.bz2
           │     │            └── 24Mg-1H__XAB.trans.bz2
           │     ├── 25Mg-1H
           │     │       ├── Yadin
           │     │       └── XAB
           │     │            ├── 25Mg-1H__XAB.def
           │     │            ├── 25Mg-1H__XAB.pf
           │     │            ├── 25Mg-1H__XAB.states.bz2
           │     │            └── 25Mg-1H__XAB.trans.bz2
           │     └── 26Mg-1H
           ├── ...
           │
```

`SavePath` is the folder path for saving all results obtained by the PyExoCross program.

*Example*

```bash
# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
```

**For HITRAN**

`ReadPath` is the file path of input line list `.par` file.

```
/home/username/data/hitran/CO2.par
```

`SavePath` is the folder path for saving all results obtained by the PyExoCross program.

*Example*

```bash
# File path #
ReadPath                                /home/jingxin/data/HITRAN/CO2.par
SavePath                                /home/jingxin/data/pyexocross/
```

## Functions

In current version, *PyExoCross* can convert data format between the ExoMol and the HITRAN formats.

*PyExoCross* also implements the computations of other useful functions including partition functions, specific heats, cooling functions, radiative lifetimes, stick spectra and cross sections for data from the ExoMol database.

*PyExoCross* provides computations of cross sections for data from the HITRAN database. If you want to use the other functions, please convert the data format from the HITRAN format to the ExoMol format at first and then treat the data as the ExoMol data to use *PyExoCross*.

Use this function or not:

`0` means no;

`1` means yes.

If the value of a function's second column is `0`, then there is no need to do any changes in this function's own section, the program won't process data with this function. Although this function won't be used by users, please don't delete this function's own section, otherwise, the program cannot run.

*Example*

```bash
# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           1
CoolingFunctions                        0
Lifetimes                               0
StickSpectra                            0
CrossSections                           1
```

## Quantum Numbers

Please provide the labels `QNslabel` and formats `QNsformat` of the quantum numbers when you use *PyExoCross* to convert data format, calculate stick spectra or cross sections if you need the quantum filter.

* The definition file `.def` of ExoMol database (available at [exomol.com](https://www.exomol.com/)) provides the labels and formats of the quantum numbers for each species for reference.
* HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.

**Note**

You can define the quantum number column name by yourself, but please make sure it has letters without any spaces.
e.g. 'c1', 'c2', 'v1', 'v2', 'electronicState', 'electronic_state', '1v', '2v', 'M/E/C'.
Wrong format of the quantum number column nams: '1', '2', 'electronic state'.

*Example*

```bash
# Quantum numbers #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %1d      %7.1f    %7.1f
```
