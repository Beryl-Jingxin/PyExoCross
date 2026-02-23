# Conversion

*PyExoCross* can convert data format between the ExoMol and HITRAN (HITRAN 2004 edition 160-character records) formats.

| Database | ExoMol format | HITRAN format |
| :------: | :-----------: | :-----------: |
| ExoMol   | &#x2705;      |               |
| ExoAtom  | &#x2705;      |               |
| HITRAN   |               | &#x2705;      |
| HITEMP   |               | &#x2705;      |

`ConversionFrequncyRange`: Provide the wavenumber range in unit of cm⁻¹.

## Data format

`ConversionFormat`

`0` means no conversion.

`1` means convert data format from ExoMol to HITRAN. In this case, `Database` should be `ExoMol`.

`2` means convert data format from HITRAN to ExoMol. In this case, `Database` should be `HITRAN`.

## Quantum number label

`GlobalQNLabel` and `LocalQNLabel`

Here, the quantum number labels are what kind of quantum numbers you want to save in the output file.

For 3 different symmetry indices and inversional parity labels, please write write them as following symbols:

`Gtot`: total symmetry index;

`Gvib`: vibrational symmetry indices;

`Grot`: rotational symmetry indices;

`taui`: inversional parity.

## Quantum number format

`GlobalQNFormat` and `LocalQNFormat`

Here, the quantum number formats are the formats of quantum numbers you want to save in the output file.

In the standard HITRAN2004 format, both global and local quantum numbers have 15 characters.

## Filters

`ConvUncFilter(Y/N)`: If it is yes, the value is the maximum uncertainty you require, in unit of cm⁻¹. 

`ConvThreshold(Y/N)`: If it is yes, the value is the minimum intensity you require, in unit of cm/molecule.

## Convert data format from ExoMol to HITRAN

**For ExoMol database**

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
SpeciesID                               501


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/MgH_ExoMol_toHITRAN.log


# Functions #
Conversion                              1
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f


# Conversion #
ConversionFormat                        1  
ConversionFrequncyRange                 0          30000          # Wavenumber in unit of cm-1    
GlobalQNLabel                           eS       v        Omega
GlobalQNFormat                          %9s      %2d      %4s
LocalQNLabel                            J        e/f
LocalQNFormat                           %5.1f    %2s
ConvUncFilter(Y/N)                      Y          0.01           # If Y, default value 0.01 cm-1
ConvThreshold(Y/N)                      Y          1e-30          # If Y, default value 1e-30 cm/molecule
```

**For ExoAtom database**

*Example*

```bash
# Data source #
Database                                ExoAtom
Atom                                    Li
Dataset                                 NIST


# File path #
ReadPath                                /mnt/data/exoatom/exoatom_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/Li_ExoAtom_toHITRAN.log


# Functions #
Conversion                              1
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               1
NCPUfiles                               1
ChunkSize                               10000


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                configuration     LS       parity  
QNsformat                               %30s              %30s     %2s      


# Conversion #
ConversionFormat                        1  
ConversionFrequncyRange                 0          43000          # Wavenumber in unit of cm-1          
GlobalQNLabel                           configuration     LS       
GlobalQNFormat                          %30s              %30s     
LocalQNLabel                            J       parity
LocalQNFormat                           %5.1f   %2s
ConvUncFilter(Y/N)                      N          0.01           # If Y, default value 0.01 cm-1
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30 cm/molecule
```

## Convert data format from HITRAN to ExoMol

**For HITRAN and HITEMP databases**

*Example*

```bash
# Data source #
Database                                HITRAN
Molecule                                NO
Isotopologue                            14N-16O
Dataset                                 XABC-HITRAN
SpeciesID                               81


# File path #
ReadPath                                /home/jingxin/data/HITRAN/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_HITRAN_toExoMol.log


# Functions #
Conversion                              1
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               16
NCPUfiles                               1
ChunkSize                               1000000


# Quantum numbers for conversion, stick spectra and cross sections 
QNslabel                                X     Omega   v1      m      Br    Sym    F
QNsformat                               %2s   %3s     %2d     %1s    %2s   %1s    %5s


# Conversion # 
ConversionFormat                        2  
ConversionFrequncyRange                 0          63000          # Wavenumber in unit of cm-1       
GlobalQNLabel                           X       Omega    v1  
GlobalQNFormat                          %2s     %3s      %2d   
LocalQNLabel                            m       Br       Sym     F   
LocalQNFormat                           %1s     %2s      %1s     %5s   
ConvUncFilter(Y/N)                      N          0.01           # If Y, default value 0.01 cm-1
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30 cm/molecule
```

**Note**

1. ExoMol format definition files `.def`, `.def.json`, and `.adef.json` (available at [exomol.com](https://www.exomol.com/)) provide the labels and formats of the quantum numbers for each species for reference.
2. HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.
