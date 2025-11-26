# Lifetimes

Please provide the line list files.

If you want to save compressed .states.bz2 file, please wirte ``Y`` after ``Compress(Y/N)``,
otherwise, if you write ``N``, the uncompressed .states file will be saved.

The lifetimes equation is:

$$ 
    \tau_i = \frac{1}{{\textstyle \sum_{f} A_{fi}}}. 
$$

*Example*

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
LogFilePath                             /home/jingxin/data/pyexocross/log/MgH_ExoMol_lifetime.log


# Functions #
Conversion                              0
PartitionFunctions                      0
SpecificHeats                           0
CoolingFunctions                        0
Lifetimes                               1
OscillatorStrengths                     0
StickSpectra                            0
CrossSections                           0


# Cores and chunks #
NCPUtrans                               4
NCPUfiles                               1 
ChunkSize                               1000000


# Calculate lifetimes #
Compress(Y/N)                           Y                         # If Y, save as .states.bz2 file; otherwise, save as .states file
```

**Note**

If the line lists data is not in the ExoMol format, please convert your
data format into the ExoMol format at first and then compute lifetime with *PyExoCross*.
So please read [**Conversion**](`https://pyexocross.readthedocs.io/en/latest/conversion.html`) and write ``1`` after ``Conversion``, ``2`` after ``ConversionFormat`` and fill ``Conversion`` section.
