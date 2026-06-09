# Partition functions

Please provide the line lists, temperature step ``Ntemp`` and the
maximum of the temperature ``Tmax``.

``Ntemp`` is always set as ``1`` K.

``Tmax`` can be set by yourself and the definition file ``.def`` from
the ExoMol database provides the maximum temperature of each molecule
for reference.

The temperatures are start from 1 K to ``Tmax`` K in the output file.

The partition functions equation is:

$$
   Q(T)=\sum_n g_n^{\textrm{tot}} e^{-c_2\tilde{E}_n/T}.
$$

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                CO2
Isotopologue                            12C-16O2
Dataset                                 UCL-4000
SpeciesID                               21


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/CO2_ExoMol_pf.log


# Functions #
Conversion                              0
PartitionFunctions                      1
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               4
ChunkSize                               1000000


# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 
```

**Note**

For HITRAN/HITEMP input, partition functions are calculated directly from the
line list records; format conversion and quantum number settings are not
required. At high temperature, the result may differ from HITRAN/TIPS values
because a `.par` file is a transition list and may not contain a complete set
of high energy states.
