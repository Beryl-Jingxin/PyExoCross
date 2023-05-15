# Examples

The examples of the whole input files for the ExoMol and HITRAN databases.

## Example for the ExoMol database

```bash
# Basic information #
Database                                ExoMol
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
mol_iso_id                              501


# File path #
ReadPath                                /mnt/data/exomol/exomol3_data/
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              0
PartitionFunctions                      0
CoolingFunctions                        0
Lifetimes                               0
SpecificHeats                           0
StickSpectra                            1
CrossSections                           0


# Quantum numbers #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %1d      %7.1f    %7.1f


# Conversion #
ConversionFormat                        1  
ConversionUncertainty                   0.01
ConversionFrequncyRange                 0                 30000  
GlobalQNLabel                           eS       v        Omega
GlobalQNFormat                          %10s     %1d      %4s
LocalQNLabel                            J        e/f
LocalQNFormat                           %5.1f    %2s
   

# Calculate partition, cooling functions or specific heats #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
None


# Calculate stick spectra or cross-sections #
Temperature                             300
Pressure                                1
Range                                   0          30000
Npoints/BinSize                         Npoints    30001

Cutoff(Y/N)                             N          25             # Default value 25
Threshold(Y/N)                          N          1e-30          # Default value 1e-30
UncFilter(Y/N)                          N          0.001
QNsFilter(Y/N)                          N          par[+]   e/f[e]   v[0,1,2,3]  
DopplerHWHM(Y/N)                        Y          3              # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     N          0.7            # Set Lorentzian HWHM as a constant

Broadeners                              Default    Air    Self    H2    He    CO2
Ratios                                  1.0        0.0    0.0     0.0   0.0   0.0

Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
Profile                                 Gaussian  
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'


Note:
1. mol_iso_id
   If the database is ExoMol: mol_iso_id   0
   If the database is HITRAN: mol_iso_id   moleculeIDisotopologueID          # e.g. 81 (NO), 261(C2H2)
2. ReadPath
   If the database is ExoMol: ReadPath  /mnt/data/exomol/exomol3_data/       # folder path
   (when file path is /mnt/data/exomol/exomol3_data/AlH/27Al-1H/AlHambra/27Al-1H__AlHambra.def)
   If the database is HITRAN: ReadPath  /home/username/data/hitran/AlH.par   # .par file path
   (when file path is /home/username/data/hitran/AlH.par)
3. Functions
   Functions part: (calculate the functions or not) 0 means no, 1 means yes. 
   Just change the information which you will use, please do not delete other information.
   (Cooling function's minimal T = 200 K, others (partition function, specific heat and lifetime) minimal T = 1 K )
4. ConversionFormat
   0: no conversion
   1: from ExoMol to HITRAN
   2: from HITRAN to ExoMol
   Note: for 3 different symmetry indices and inversional parity labels, please write write them as following symbols:
         Gtot: total symmetry index;
         Gvib: vibrational symmetry indices;
         Grot: rotational symmetry indices;
         taui: inversional parity.
5. Broadeners; Ratios
   The broadening types and ratio values are corresponding, please do not change the place or delete elements. 
6. Range
   Give two values as the minimum and maximum of the wavenumber range. No ',' or ';' between these two numbers, just leave spaces here.
7. Npoints/BinSize
   e.g. 'Npoints      100001' or 'BinSize    1'.
8. Cutoff(Y/N); Threshold(Y/N); UncFilter(Y/N); QNsFilter(Y/N)
   'Y', 'YES', 'Yes', 'yes' and 'N', 'NO', 'No', 'no' all can work.
   (1). Cutoff(Y/N)   : e.g. 'Y          25'    or 'N';
   (2). Threshold(Y/N): e.g. 'Y          1e-30' or 'N';
   (3). UncFilter(Y/N): e.g. 'Y          0.001' or 'N';
   (4). QNsFilter(Y/N): e.g. 'Y          +/-[+]   e/f[e]   v[0,1,2,3]  ' or 'N';
9. Quantum number filter
   (1). QNslabel                         +/-  e/f   eS    v     Lambda   Sigma    Omega
   (2). QNsformat                        %1s  %1s   %13s  %3d   %1d      %7.1f    %7.1f
   Note: you can define the QN column name by yourself, but please make sure it has letters without any spaces.
   e.g. 'c1'  'c2'  'v1'  'v2'  'electronicState'  'electronic_state'  '1v'  '2v'  'M/E/C'.
   Wrong format of the QN column nams: '1'   '2'   'electronic state'.
10. Profile
   Choose line profile from: 
   Doppler, Gaussian, Lorentzian, SciPyVoigt, SciPyWofzVoigt, 
   PseudoVoigt, PseudoKielkopfVoigt, PseudoOliveroVoigt, PseudoLiuLinVoigt, PseudoRoccoVoigt,
   BinnedDoppler, BinnedGaussian, BinnedLorentzian, BinnedVoigt.
```

## Example for the HITRAN database

```bash
# Basic information #
Database                                HITRAN
Molecule                                MgH
Isotopologue                            24Mg-1H
Dataset                                 XAB
mol_iso_id                              501


# File path #
ReadPath                                /home/jingxin/data/pyexocross/conversion/24Mg-1H__XAB.par
SavePath                                /home/jingxin/data/pyexocross/


# Functions #
Conversion                              1
PartitionFunctions                      0
CoolingFunctions                        0
Lifetimes                               0
SpecificHeats                           0
StickSpectra                            0
CrossSections                           0


# Quantum numbers #
QNslabel                                par  e/f   eS    v     Lambda   Sigma    Omega
QNsformat                               %1s  %1s   %13s  %3d   %1d      %7.1f    %7.1f


# Conversion #
ConversionFormat                        2  
ConversionUncertainty                   0.005
ConversionFrequncyRange                 0                 30000  
GlobalQNLabel                           eS       v        Omega
GlobalQNFormat                          %10s     %1d      %4s
LocalQNLabel                            J        e/f
LocalQNFormat                           %5.1f    %2s
   

# Calculate partition, cooling functions or specific heats #
Ntemp                                   1                         # The number of temperature steps
Tmax                                    5000                      # Maximal temperature in K 


# Calculate lifetimes #
None


# Calculate stick spectra or cross-sections #
Temperature                             300
Pressure                                1
Range                                   0          30000
Npoints/BinSize                         Npoints    30001

Cutoff(Y/N)                             N          25             # Default value 25
Threshold(Y/N)                          N          1e-30          # Default value 1e-30
UncFilter(Y/N)                          N          0.001
QNsFilter(Y/N)                          N          par[+]   e/f[e]   v[0,1,2,3]  
DopplerHWHM(Y/N)                        Y          3              # Set Doppler HWHM as a constant
LorentzianHWHM(Y/N)                     Y          0.7            # Set Lorentzian HWHM as a constant

Broadeners                              Default    Air    Self    H2    He    CO2
Ratios                                  1.0        0.0    0.0     0.0   0.0   0.0

Absorption/Emission                     Absorption                # 'Absorption' or 'Emission'
Profile                                 SciPyVoigt  
Wavenumber(wn)/wavelength(wl)           wn                        # 'wn' or 'wl'


Note:
1. mol_iso_id
   If the database is ExoMol: mol_iso_id   0
   If the database is HITRAN: mol_iso_id   moleculeIDisotopologueID          # e.g. 81 (NO), 261(C2H2)
2. ReadPath
   If the database is ExoMol: ReadPath  /mnt/data/exomol/exomol3_data/       # folder path
   (when file path is /mnt/data/exomol/exomol3_data/AlH/27Al-1H/AlHambra/27Al-1H__AlHambra.def)
   If the database is HITRAN: ReadPath  /home/username/data/hitran/AlH.par   # .par file path
   (when file path is /home/username/data/hitran/AlH.par)
3. Functions
   Functions part: (calculate the functions or not) 0 means no, 1 means yes. 
   Just change the information which you will use, please do not delete other information.
   (Cooling function's minimal T = 200 K, others (partition function, specific heat and lifetime) minimal T = 1 K )
4. ConversionFormat
   0: no conversion
   1: from ExoMol to HITRAN
   2: from HITRAN to ExoMol
   Note: for 3 different symmetry indices and inversional parity labels, please write write them as following symbols:
         Gtot: total symmetry index;
         Gvib: vibrational symmetry indices;
         Grot: rotational symmetry indices;
         taui: inversional parity.
5. Broadeners; Ratios
   The broadening types and ratio values are corresponding, please do not change the place or delete elements. 
6. Range
   Give two values as the minimum and maximum of the wavenumber range. No ',' or ';' between these two numbers, just leave spaces here.
7. Npoints/BinSize
   e.g. 'Npoints      100001' or 'BinSize    1'.
8. Cutoff(Y/N); Threshold(Y/N); UncFilter(Y/N); QNsFilter(Y/N)
   'Y', 'YES', 'Yes', 'yes' and 'N', 'NO', 'No', 'no' all can work.
   (1). Cutoff(Y/N)   : e.g. 'Y          25'    or 'N';
   (2). Threshold(Y/N): e.g. 'Y          1e-30' or 'N';
   (3). UncFilter(Y/N): e.g. 'Y          0.001' or 'N';
   (4). QNsFilter(Y/N): e.g. 'Y          +/-[+]   e/f[e]   v[0,1,2,3]  ' or 'N';
9. Quantum number filter
   (1). QNslabel                         +/-  e/f   eS    v     Lambda   Sigma    Omega
   (2). QNsformat                        %1s  %1s   %13s  %3d   %1d      %7.1f    %7.1f
   Note: you can define the QN column name by yourself, but please make sure it has letters without any spaces.
   e.g. 'c1'  'c2'  'v1'  'v2'  'electronicState'  'electronic_state'  '1v'  '2v'  'M/E/C'.
   Wrong format of the QN column nams: '1'   '2'   'electronic state'.
10. Profile
   Choose line profile from: 
   Doppler, Gaussian, Lorentzian, SciPyVoigt, SciPyWofzVoigt, 
   PseudoVoigt, PseudoKielkopfVoigt, PseudoOliveroVoigt, PseudoLiuLinVoigt, PseudoRoccoVoigt,
   BinnedDoppler, BinnedGaussian, BinnedLorentzian, BinnedVoigt.
```
