# Specific heats

Please provide the line lists, temperature step ``Ntemp`` and the
maximum of the temperature ``Tmax``.

``Ntemp`` is always set as ``1`` K.

``Tmax`` can be set by yourself and the definition file ``.def`` from
the ExoMol database provides the maximum temperature of each molecule
for reference.

The temperatures are start from 1 K to ``Tmax`` K in the output file.

The specific heats equation is:

$$
   C_p(T) = R\left [\frac{Q''}{Q}-\left (\frac{Q'}{Q} \right )^2 \right ]+\frac{5R}{2},
$$

where the partition function :math:`Q(T)` and its first two moments are:

$$
   Q(T)=\sum_n g_n^{\textrm{tot}}e^{-c_2\tilde{E}_n/T}, 
$$

$$
   Q'(T) = T\frac{\mathrm{d} Q}{\mathrm{d} T} =\sum_n 
   g_n^{\textrm{tot}}\left(\frac{c_2 \tilde{E}_n}{T}\right)\exp\left(-\frac{c_2 \tilde{E}_n}{T}\right),
$$

$$
   Q''(T) = T^2\frac{\mathrm{d}^2 Q}{\mathrm{d} T^2}+2Q' =\sum_n g_n^{\textrm{tot}}
   \left(\frac{c_2 \tilde{E}_n}{T}\right)^2\exp\left(\frac{c_2 \tilde{E}_n}{T}\right).
$$

*Example: ExoMol*

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
LogFilePath                             /home/jingxin/data/pyexocross/log/CO2_ExoMol_cp.log


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           1
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               1
NCPUfiles                               1
ChunkSize                               1000000


# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 
```

*Example: HITRAN*

```bash
# Data source #
Database                                HITRAN
Molecule                                NO
Isotopologue                            14N-16O
Dataset                                 NO-HITRAN
SpeciesID                               81


# File path #
ReadPath                                /home/jingxin/data/HITRAN/
SavePath                                /home/jingxin/data/pyexocross/
LogFilePath                             /home/jingxin/data/pyexocross/log/NO_HITRAN_cp.log


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           1
CoolingFunctions                        0
Lifetimes                               0
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               1
NCPUfiles                               1
ChunkSize                               1000000


# Calculate partition, specific heats or cooling functions #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K
```
