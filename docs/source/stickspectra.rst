Stick spectra
=============

Please provide the line lists, ``Temperature``, wavenumber ``Range``, 
uncertainty filter ``UncFilter`` and ``PlotStickSpectra(Y/N)``.

Filters
:::::::

``Y/N``: ``Y``, ``YES``, ``Yes``, ``yes`` and ``N``, ``NO``, ``No``, ``no`` all can work. \
If you don't use it, write ``N`` here. You don't need to change the content behind it. \
If after using filters, the program gets an empty result, then you will receive an warning to ask you write new filter values.

If ``UncFilter(Y/N)`` is yes, the value is the maximum uncertainty you require. 

If ``Threshold(Y/N)`` is yes, the value is the minimum intensity you require.

If ``QNsFilter(Y/N)`` is yes, program will do filter on quantum numbers. \
Write the quantum number labels required here, the spelling of the quantum number labels must be the same as ``QNslabel``. \
The other quantum number labels which are in ``QNslabel`` but not in ``QNsFilter(Y/N)`` will not be stored in the result file. \
Write the quantum number values required after each label in ``[]`` and seperated by ``,``, don't leave blank between different values inside the ``[]``. \
Leave blank between different quantum number labels, don't write ``,``.\
If you need all values of a quantum number label, write this label and wirte nothing inside the ``[]``. Note, don't write any blank inside ``[]``, you should write ``[]``, not ``[ ]``.\
Inside ``[]`` use ``,``, don't write any blank inside the ``[]``. Outside ``[]``, use blank, don't write any ``,`` outside ``[]``.

If you need a stick spectra figure, please write ``Y`` after ``PlotStickSpectra(Y/N)``.

The intensity equation is:

.. math::

    I(f \gets i) = \frac{g'{A}_{fi}}{8 \pi c \tilde{v}^2_{fi}} 
    \frac{e^{-c_2 \tilde{E}'' / T} (1 - e^{-c_2 \tilde{v}_{fi} 
    / T })}{Q(T)}.

*Example*

.. code:: bash

    # Quantum numbers for conversion, stick spectra and cross sections #
    QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
    QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f

    # Calculate stick spectra or cross sections #
    Temperature                             300
    Range                                   0          30000
    Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
    UncFilter(Y/N)                          Y          0.001          # If Y, default value 0.001
    Threshold(Y/N)                          Y          1e-30          # If Y, default value 1e-30
    QNsFilter(Y/N)                          Y          par[]   e/f[e]   v[0,1,2,3]  

    # Calculate stick spectra #
    PlotStickSpectra(Y/N)                   Y

.. code:: bash

    # Quantum numbers for conversion, stick spectra and cross sections #
    QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
    QNsformat                               %1s  %1s   %13s  %3d   %2d      %7.1f    %7.1f

    # Calculate stick spectra or cross sections #
    Temperature                             1000
    Range                                   1000       5000
    Absorption/Emission                     Emission                  # 'Absorption' or 'Emission'
    UncFilter(Y/N)                          N          0.005          # If Y, default value 0.001
    Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
    QNsFilter(Y/N)                          N          par[]   e/f[e]   v[0,1,2,3]  

    # Calculate stick spectra #
    PlotStickSpectra(Y/N)                   N
