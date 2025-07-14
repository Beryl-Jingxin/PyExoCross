# Filters

`Y/N`: `Y`, `YES`, `Yes`, `yes` and `N`, `NO`, `No`, `no` all can work. \
If you don't use it, write `N` here. You don't need to change the content behind it. \
If after using filters, the program gets an empty result, then you will receive an warning to ask you write new filter values to icrease the range.

## Uncertainty filters

If `UncFilter(Y/N)` is yes, the value is the maximum uncertainty you require, in unit $\textrm{cm}^{-1}$. 

## Threshold filters

If `Threshold(Y/N)` is yes, the value is the minimum intensity you require, in unit cm/molecule.

## Quantum number filters

If `QNsFilter(Y/N)` is yes, program will do filter on quantum numbers. 

### Quantum number label filters

Write the quantum number labels required here, the spelling of the quantum number labels must be the same as `QNslabel`. 
The other quantum number labels which are in `QNslabel` but not in `QNsFilter(Y/N)` will not be stored in the result file. 

### Quantum number value filters

Write the quantum number values required after each label in `[]` and seperated by `,` and `;`, don't leave blank between different values inside the `[]`. \
Leave blank between different quantum number labels, don't write `,`. \
If you need all values of a quantum number label, write this label and wirte nothing inside the `[]`. Note, don't write any blank inside `[]`, you should write `[]`, not `[ ]`. \
Inside `[]` use `,` and `;` don't write any blank inside the `[]`. Outside `[]`, use blank , don't write any `,` or `;` outside `[]`. \
For one quantum number label, write in one `[]`, you can provide the quantum number values for upper and lower states, and seperated by `,`. \
For one quantum number label, write in one `[]`, you can provide more than one pair of values, and seperated by `;`.

*Example*

`v1[]` means you want quantum number label v1 and you want all quantum number values of this label v1. \
`v1[1,0]` means you want quantum number label v1 and you want the upper QN = 1 and lower QN = 0. So v1' = 1 and v1" = 0. \
`v1[,0]` means you want quantum number label v1 and you want all upper QN but the lower QN = 0. So v1' = 0, 1, 2, 3, 4, ... and v1" = 0. \
`v1[3,]` means you want quantum number label v1 and you want all lower QN but the upper QN = 3. So v1' = 3 and v1" = 0, 1, 2, 3, 4, ... \
`v1[1,1;2,2]` means you want quantum number label v1 and you want when v1' = 1, v1" = 1; when v1' = 2, v1" = 2. \
`v1[1,;,0;5,5]  v2[]` means you want quantum number labels v1 and v2. For v1, you want all lines with v1' = 1 , all lines with v1" = 0 and the lines with v1' = 5 and at the same time v1" = 5. Meanwhile, you want all lines for v2.

* The definition file `.def`, `.def.json`, or `.adef.json` of ExoMol database format (available at [exomol.com](https://www.exomol.com/)) provides the labels and formats of the quantum numbers for each species for reference.
* HITRAN2020 supplementary material ([link](https://hitran.org/media/refs/HITRAN_QN_formats.pdf)) provides the notation and format for quanta identifications for reference.

**Note**

You can define the quantum number column name by yourself, but please make sure it has letters without any blanks.
e.g. 'c1', 'c2', 'v1', 'v2', 'electronicState', 'electronic_state', '1v', '2v', 'M/E/C'.
Wrong format of the quantum number column nams: '1', '2', 'electronic state'.

*Example*

```bash
UncFilter(Y/N)                          N          0.01           # If Y, default value 0.01
Threshold(Y/N)                          N          1e-30          # If Y, default value 1e-30
QNsFilter(Y/N)                          N          par[]   e/f[]   v[1,;2,2;2,1;,0]  
```
