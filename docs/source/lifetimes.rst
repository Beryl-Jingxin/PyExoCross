Lifetimes
===================

| Please provide the line list files. 
| If you want to save compressed .states.bz2 file, please wirte ``Y`` after ``Compress(Y/N)``, 
otherwise, if you write ``N``, the uncompressed .states file will be saved.

The lifetimes equation is:

.. math::

   \tau_i = \frac{1}{{\textstyle \sum_{f} A_{fi}}}.

*Example*

.. code:: bash

    # Calculate lifetimes #
    Compress(Y/N)                           Y                         If Y, save as .states.bz2 file; otherwise, save as .states file


**Note**

If the line lists data is not in the ExoMol format, please convert your
data format into the ExoMol format at first and then compute lifetime with *PyExoCross*.
So please read `**Conversion** <https://pyexocross.readthedocs.io/en/latest/conversion.html>`_ 
and write ``1`` after ``Conversion``, ``2`` after ``ConversionFormat`` and fill ``Conversion`` section.
 