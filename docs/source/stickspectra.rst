Stick spectra
=============

Please provide the line lists, ``Temperature``, ``Pressure``, 
wavenumber ``Range`` and uncertainty filter ``UncFilter``.

If you don't use the uncertainty filter ``UncFilter``, write ``N`` here. 
You don't need to change the number behind it.

The intensity equation is:

.. math::

    I(f \gets i) = \frac{g'{A}_{fi}}{8 \pi c \tilde{v}^2_{fi}} 
    \frac{e^{-c_2 \tilde{E}'' / T} (1 - e^{-c_2 \tilde{v}_{fi} 
    / T })}{Q(T)}.

*Example*

.. code:: bash

    # Calculate stick spectra or cross-sections #
    Temperature                             300
    Pressure                                1
    Range                                   0          30000

    UncFilter(Y/N)                          Y          0.001

.. code:: bash
    
    # Calculate stick spectra or cross-sections #
    Temperature                             1000
    Pressure                                1.5
    Range                                   1000       5000

    UncFilter(Y/N)                          N          0.005


**Note**

If the line lists data is not in the ExoMol format, 
please convert your data format into the ExoMol format at first 
and then compute partition functions with *PyExoCross*.
