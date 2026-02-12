# Cooling functions

Please provide the line lists, temperature step ``Ntemp``
and the maximum of the temperature ``Tmax``.

``Ntemp`` is always set as ``1`` K.

``Tmax`` can be set by yourself and the definition file ``.def`` from
the ExoMol database provides the maximum temperature of each molecule
for reference.

The temperatures are start from 1 K to ``Tmax`` K in the output file.

The cooling functions equation is:

$$ 
   W(T) = \frac{1}{4 \pi Q(T)} \sum_{f,i} A_{fi} h c \tilde{v}_{fi} g' e^{-c_2 \tilde{E}' / T}, 
$$

where $Q(T)$ is the partition function.

$$ 
   Q(T)=\sum_n g_n^{\textrm{tot}} e^{-c_2\tilde{E}_n/T}. 
$$

*Example*

```bash
# Data source #
Database                                ExoMol
Molecule                                H2O
Isotopologue                            1H2-16O
Dataset                                 POKAZATEL
SpeciesID                               11


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/H2O_ExoMol_cf.log


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        1
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               32
ChunkSize                               100000


# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    10000                     # Maximal temperature in K 
```
