Stick spectra
=============

Please provide the line lists, ``Temperature``, wavenumber ``Range``, 
uncertainty filter ``UncFilter`` and ``PlotStickSpectra(Y/N)``.

If you need a stick spectra figure, please write ``Y`` after ``PlotStickSpectra(Y/N)``.
If you don't use the uncertainty filter ``UncFilter``, write ``N`` here. 
You don't need to change the number behind it.

The intensity equation is:

.. math::

    I(f \gets i) = \frac{g'{A}_{fi}}{8 \pi c \tilde{v}^2_{fi}} 
    \frac{e^{-c_2 \tilde{E}'' / T} (1 - e^{-c_2 \tilde{v}_{fi} 
    / T })}{Q(T)}.

*Example*

.. code:: bash

    # Calculate stick spectra or cross sections #
    Temperature                             300
    Range                                   0          30000
    Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
    UncFilter(Y/N)                          Y          0.001          # If Y, default value 0.001
    Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30

    # Calculate stick spectra #
    PlotStickSpectra(Y/N)                   Y

.. code:: bash
    
    # Calculate stick spectra or cross sections #
    Temperature                             1000
    Range                                   1000       5000
    Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
    UncFilter(Y/N)                          N          0.005          # If Y, default value 0.001
    Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30

    # Calculate stick spectra #
    PlotStickSpectra(Y/N)                   N

**Note**

If the line lists data is not in the ExoMol format, 
please convert your data format into the ExoMol format at first 
and then compute partition functions with *PyExoCross*.
