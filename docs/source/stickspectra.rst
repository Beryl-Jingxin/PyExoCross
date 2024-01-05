Stick spectra
=============

``Temperature``: Please provide temperature in unit K.

``Range``: Give two values as the minimum and maximum of the wavenumber range. No ``,`` or ``;`` 
between these two numbers, just leave blank here.

``Absorption/Emission``: Choose ``Absorption`` or ``Emission``.

If you want a figure of corss sections, please set ``Y`` for ``PlotCrossSection(Y/N)``.

And if you want set the lower limit of y-axis for plotting, please write after ``Y-axisLimitStick``, 
otherwise, the default lower limit y-axis is 1e-30.

Filters
:::::::

| ``Y/N``: ``Y``, ``YES``, ``Yes``, ``yes`` and ``N``, ``NO``, ``No``, ``no`` all can work. 
| If you don't use it, write ``N`` here. You don't need to change the content behind it. 
| If after using filters, the program gets an empty result, then you will receive an warning to ask you write new filter values to increase the range. 

If ``UncFilter(Y/N)`` is yes, the value is the maximum uncertainty you require. 

If ``Threshold(Y/N)`` is yes, the value is the minimum intensity you require.

| If ``QNsFilter(Y/N)`` is yes, program will do filter on quantum numbers.
| Write the quantum number labels required here, the spelling of the quantum number labels must be the same as ``QNslabel``. 
| The other quantum number labels which are in ``QNslabel`` but not in ``QNsFilter(Y/N)`` will not be stored in the result file. 
| Write the quantum number values required after each label in ``[]`` and seperated by ``,`` and ``;``, don't leave blank between different values inside the ``[]``. 
| Leave blank between different quantum number labels, don't write ``,``.
| If you need all values of a quantum number label, write this label and wirte nothing inside the ``[]``. Note, don't write any blank inside ``[]``, you should write ``[]``, not ``[ ]``.
| Inside ``[]`` use ``,`` and ``;``, don't write any blank inside the ``[]``. Outside ``[]``, use blank, don't write any ``,`` or ``;`` outside ``[]``. 
| For one quantum number label, write in one ``[]``, you can provide the quantum number values for upper and lower states, and seperated by ``,``. 
| For one quantum number label, write in one ``[]``, you can provide more than one pair of values, and seperated by ``;``.

*Example*

| ``v1[]`` means you want quantum number label v1 and you want all quantum number values of this label v1.
| ``v1[1,0]`` means you want quantum number label v1 and you want the upper QN = 1 and lower QN = 0. So v1' = 1 and v1" = 0.
| ``v1[,0]`` means you want quantum number label v1 and you want all upper QN but the lower QN = 0. So v1' = 0, 1, 2, 3, 4, ... and v1" = 0. 
| ``v1[3,]`` means you want quantum number label v1 and you want all lower QN but the upper QN = 3. So v1' = 3 and v1" = 0, 1, 2, 3, 4, ... 
| ``v1[1,1;2,2]`` means you want quantum number label v1 and you want when v1' = 1, v1" = 1; when v1' = 2, v1" = 2.
| ``v1[1,;,0;5,5]  v2[]`` means you want quantumnumber labels v1 and v2. For v1, you want all lines with v1' = 1 , all lines with v1" = 0 and the lines with v1' = 5 and at the same time v1" = 5. Meanwhile, you want all lines for v2.


If you need a stick spectra figure, please write ``Y`` after ``PlotStickSpectra(Y/N)``. 

The intensity equation is:

.. math::

    I(f \gets i) = \frac{g'{A}_{fi}}{8 \pi c \tilde{v}^2_{fi}} 
    \frac{e^{-c_2 \tilde{E}'' / T} (1 - e^{-c_2 \tilde{v}_{fi} 
    / T })}{Q(T)}.

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
    OscillatorStrengths                     0
    StickSpectra                            1
    CrossSections                           0


    # Quantum numbers for conversion, stick spectra and cross sections #
    QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
    QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f


    # Calculate stick spectra or cross sections #
    Temperature                             300
    Range                                   0          30000
    Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
    UncFilter(Y/N)                          Y          0.001          # If Y, default value 0.01
    Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30
    QNsFilter(Y/N)                          Y          par[]   e/f[]   v[1,1;1,0;2,;,0]  


    # Calculate stick spectra #
    PlotStickSpectra(Y/N)                   Y
    Y-axisLimitStick                        1e-40                     # Default value is 1e-30

.. code:: bash

    # Data source #
    Database                                HITRAN
    Molecule                                NO
    Isotopologue                            14N-16O
    Dataset                                 NO-HITRAN
    MolIsoID                                81


    # File path #
    ReadPath                                /home/jingxin/data/HITRAN/NO.par
    SavePath                                /home/jingxin/data/pyexocross/


    # Functions #
    Conversion                              0
    PartitionFunctions                      0
    SpecificHeats                           0
    CoolingFunctions                        0
    Lifetimes                               0
    OscillatorStrengths                     0
    StickSpectra                            1
    CrossSections                           0


    # Quantum numbers for conversion, stick spectra and cross sections #
    QNslabel                                J       X     Omega   v1      m      Sym    
    QNsformat                               %5s     %2s   %3s     %2d     %1s    %1s


    # Calculate stick spectra or cross sections #
    Temperature                             1000
    Range                                   1000       5000
    Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
    UncFilter(Y/N)                          N          0.005          # If Y, default value 0.01
    Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
    QNsFilter(Y/N)                          N          par[]   e/f[e,e]   v[1,1;1,0]  


    # Calculate stick spectra #
    PlotStickSpectra(Y/N)                   N
    Y-axisLimitStick                                                  # Default value is 1e-30
