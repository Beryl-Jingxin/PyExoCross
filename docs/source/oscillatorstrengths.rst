Oscillator strengths
====================

Oscillator strength ``fg``, which is what people usually use.  
Sometimes people need ``f``, not ``fg``. ``f`` is the actual oscillator strength and the actual value of ``f`` is ``fg`` divided by ``g"``.

| ``fg/f``: Choose oscillator strength ``fg`` or actual oscillator strength ``f``.
| ``Ncolumns``: Type ``3`` or ``4``. 
| ``3`` columns: upper ID, lower ID, oscillator strength (ExoMol database) or g', g", oscillator strength (HITRAN database).
| ``4`` columns: upper ID, lower ID, Einstein A-coefficient, oscillator strength (ExoMol database) or g', g", Einstein A-coefficient, oscillator strength (HITRAN database).

The oscillator strengths equation is:

.. math::

    fg^{''}=\frac{g^{'}_\textrm{tot}A_{fi}}{(c\tilde{v}_{fi})^2}.

*Example*

.. code:: bash

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
    Conversion                              0
    PartitionFunctions                      0
    SpecificHeats                           0
    CoolingFunctions                        0
    Lifetimes                               0
    OscillatorStrengths                     1
    StickSpectra                            0
    CrossSections                           0


    # Calculate oscillator strengths #
    fg/f                                    fg
    Ncolumns                                4                         # 3 (without A) or 4 (with A)

