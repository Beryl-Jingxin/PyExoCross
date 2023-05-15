Specific heats
==============

Please provide the line lists, temperature step ``Ntemp`` and the
maximum of the temperature ``Tmax``.

``Ntemp`` is always set as ``1`` K.

``Tmax`` can be set by yourself and the definition file ``.def`` from
the ExoMol database provides the maximum temperature of each molecule
for reference.

The temperatures are start from 1 K to ``Tmax`` K in the output file.

The specific heats equation is:

.. math::

   C_p(T) = R\left [\frac{Q''}{Q}-\left (\frac{Q'}{Q} \right )^2 \right ]+\frac{5R}{2},

where the partition function :math:`Q(T)` and its first two moments are:

.. math::

   Q(T)=\sum_n g_n^{\textrm{tot}}e^{-c_2\tilde{E}_n/T}, 

.. math::
   
   Q’(T) = T\frac{\mathrm{d} Q}{\mathrm{d} T} =\sum_n 
   g_n^{\textrm{tot}}\left(\frac{c_2 \tilde{E}_n}{T}\right)\exp\left(-\frac{c_2 \tilde{E}_n}{T}\right),

.. math::
   
   Q’‘(T) = T^2\frac{\mathrm{d}^2 Q}{\mathrm{d} T^2}+2Q’ =\sum_n g_n^{\texttt{tot}}
   \left(\frac{c_2 \tilde{E}_n}{T}\right)^2\exp\left(\frac{c_2 \tilde{E}_n}{T}\right).

*Example*

.. code:: bash
   
   # Calculate partition, specific heats or cooling functions #
   Ntemp                                   1                         # The number of temperature steps
   Tmax                                    5000                      # Maximal temperature in K 

**Note**

If the line lists data is not in the ExoMol format, please convert your
data format into the ExoMol format at first and then compute specific
heats with *PyExoCross*.
