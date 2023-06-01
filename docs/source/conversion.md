# Conversion

*PyExoCross* can convert data format between the ExoMol and HITRAN formats.

## Data format

`ConversionFormat`

`0` means no conversion.

`1` means convert data format from ExoMol to HITRAN. In this case, `Database` should be `ExoMol`.

`2` means convert data format from HITRAN to ExoMol. In this case, `Database` should be `HITRAN`.

## Quantum number label

`GlobalQNLabel` and `LocalQNLabel`

Here, the quantum number labels are what kind of quantum numbers you want to save in the output file.

For 3 different symmetry indices and inversional parity labels, please write write them as following symbols:

`Gtot`: total symmetry index;

`Gvib`: vibrational symmetry indices;

`Grot`: rotational symmetry indices;

`taui`: inversional parity.

## Quantum number format

`GlobalQNFormat` and `LocalQNFormat`

Here, the quantum number formats are the formats of quantum numbers you want to save in the output file.

In the standard HITRAN format, both global and local quantum numbers have 15 characters.

*Example*

```bash
# Conversion #
ConversionFormat                        1  
ConversionUncertainty                   0.01
ConversionFrequncyRange                 0                 30000  
GlobalQNLabel                           eS       v        Omega
GlobalQNFormat                          %10s     %1d      %4s
LocalQNLabel                            J        e/f
LocalQNFormat                           %5.1f    %2s
```

**Note**

1. ExoMol definition file `.def` (available at [exomol.com](https://www.exomol.com/)) provides the labels and formats of the quantum numbers for each species for reference.
2. HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.
