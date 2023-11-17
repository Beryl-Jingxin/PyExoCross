# Conversion

*PyExoCross* can convert data format between the ExoMol and HITRAN formats.

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

In the standard HITRAN format, both global and local quantum numbers have 15 characters.

*Example*

Convert data format from ExoMol to HITRAN.

```bash
# Data source #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
MolIsoID                                501


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              1
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Quantum numbers for conversion, stick spectra and cross sections #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f


# Conversion #
ConversionFormat                        1  
ConversionFrequncyRange                 0                 30000  
GlobalQNLabel                           eS       v        Omega
GlobalQNFormat                          %9s      %2d      %4s
LocalQNLabel                            J        e/f
LocalQNFormat                           %5.1f    %2s
ConvUncFilter(Y/N)                      Y          0.01           # If Y, default value 0.01
ConvThreshold(Y/N)                      Y          1e-30          # If Y, default value 1e-30
```

Convert data from HITRAN to ExoMol.

```bash
# Data source #
Database                                HITRAN
Molecule                                NO
Isotopologue                            14N-16O
Dataset                                 XABC-HITRAN
MolIsoID                                81


# File path #
ReadPath                                /home/jingxin/data/HITRAN/NO.par
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              1
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Quantum numbers for conversion, stick spectra and cross sections 
QNslabel                                X     Omega   v1      m      Br    Sym    F
QNsformat                               %2s   %3s     %2d     %1s    %2s   %1s    %5s


# Conversion # 
ConversionFormat                        2  
ConversionFrequncyRange                 0          63000   
GlobalQNLabel                           X       Omega    v1  
GlobalQNFormat                          %2s     %3s      %2d   
LocalQNLabel                            m       Br       Sym     F   
LocalQNFormat                           %1s     %2s      %1s     %5s   
ConvUncFilter(Y/N)                      N          0.01           # If Y, default value 0.01
ConvThreshold(Y/N)                      N          1e-30          # If Y, default value 1e-30  
```

**Note**

1. ExoMol definition file `.def` (available at [exomol.com](https://www.exomol.com/)) provides the labels and formats of the quantum numbers for each species for reference.
2. HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.
