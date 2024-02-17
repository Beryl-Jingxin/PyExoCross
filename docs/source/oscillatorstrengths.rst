Oscillator strengths
====================

Weighted oscillator strength ``gf``, which is what people usually use.  
``f`` is the actual oscillator strength and the actual value of ``f`` is ``gf`` divided by ``g"``.
Sometimes people need oscillator strength ``f``, not ``gf``. 

The ressult file have 4 columns: 
ExoMol data: upper id, lower id, oscillator strength and wavenumber v.
HITRAN data: g', g", oscillator strength and wavenuber v.

| ``gf/f``: Choose weigthed oscillator strength ``gf`` or actual oscillator strength ``f``.
| ``PlotOscillatorStrength(Y/N)``: If you need a oscillator strength figure, please write ``Y`` here. 
| ``Y-axisLimitOscillatorStrength``: If you want set the lower limit of y-axis for plotting, please write here, otherwise, the default lower limit y-axis is 1e-30.

The oscillator strengths equation is:

.. math::

    gf=\frac{g^{'}_\textrm{tot}A_{fi}}{(c\tilde{v}_{fi})^2}.

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
    PartitionFunctions                      0
    SpecificHeats                           0
    CoolingFunctions                        0
    Lifetimes                               0
    OscillatorStrengths                     1
    StickSpectra                            0
    CrossSections                           0


    # Cores and chunks #
    NCPU                                    32
    ChunkSize                               1000000


    # Calculate oscillator strengths #
    gf/f                                    gf
    PlotOscillatorStrength(Y/N)             N    
    Y-axisLimitOscillatorStrength           1e-30                     # Default value is 1e-30

