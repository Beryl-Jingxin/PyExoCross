Partition functions
===================

Please provide the line lists, temperature step ``Ntemp`` and the
maximum of the temperature ``Tmax``.

``Ntemp`` is always set as ``1`` K.

``Tmax`` can be set by yourself and the definition file ``.def`` from
the ExoMol database provides the maximum temperature of each molecule
for reference.

The temperatures are start from 1 K to ``Tmax`` K in the output file.

The partition functions equation is:

.. math::

   Q(T)=\sum_n g_n^{\textrm{tot}} e^{-c_2\tilde{E}_n/T}.

*Example*

.. code:: bash

   # Data source #
   Database                                ExoMol
   Molecule                                CO2
   Isotopologue                            12C-16O2
   Dataset                                 UCL-4000
   MolIsoID                                21


   # File path #
   ReadPath                                /mnt/data/exomol/exomol3_data/
   SavePath                                /home/jingxin/data/pyexocross/


   # Functions #
   Conversion                              0
   PartitionFunctions                      1
   SpecificHeats                           0
   CoolingFunctions                        0
   Lifetimes                               0
   OscillatorStrengths                     0
   StickSpectra                            0
   Non-LTE                                 0
   CrossSections                           0


   # Cores and chunks #
   NCPUtrans                               32
   NCPUfiles                               32
   ChunkSize                               1000000


   # Calculate partition, specific heats or cooling functions #
   Ntemp                                   1                         # The number of temperature steps
   Tmax                                    5000                      # Maximal temperature in K 

**Note**

If the line lists data is not in the ExoMol format, please convert your
data format into the ExoMol format at first and then compute partition
functions with *PyExoCross*. 
So please read `**Conversion** <https://pyexocross.readthedocs.io/en/latest/conversion.html>`_ 
and write ``1`` after ``Conversion``, ``2`` after ``ConversionFormat`` and fill ``Conversion`` section.
