# Error report

There are some issues to solve the problems which you may meet when the program gives the following possible error reports.

###### Please type the correct uncertainty filter choice 'Y' or 'N' into the input file.

You need to enter `Y`, `YES`, `Yes`, `y`, `yes` or `N`, `NO`, `No`, `n`, `no` in the input file.

###### Please type the correct threshold choice 'Y' or 'N' into the input file.

You need to enter `Y`, `YES`, `Yes`, `y`, `yes` or `N`, `NO`, `No`, `n`, `no` in the input file.

###### Please type the correct quantum number filter choice 'Y' or 'N' into the input file.

You need to enter `Y`, `YES`, `Yes`, `y`, `yes` or `N`, `NO`, `No`, `n`, `no` in the input file.

###### Please type the correct grid choice 'Npoints' or 'BinSize' into the input file.

You need to enter `Npoints` or `BinSize`  and corresponding value in the input file.

###### Please type the correct predissociative choice 'Y' or 'N' into the input file.

You need to enter `Y`, `YES`, `Yes`, `y`, `yes` or `N`, `NO`, `No`, `n`, `no` in the input file.

###### Please type the correct cutoff choice 'Y' or 'N' into the input file.

You need to enter `Y`, `YES`, `Yes`, `y`, `yes` or `N`, `NO`, `No`, `n`, `no` in the input file.

###### Gaussian line profile requires a HWHM. Please choose 'Y' and give a value for Doppler HWHM in the input file. Otherwise, please choose Doppler line profile (with calculated temperature-dependent Doppler HWHM).

Gaussian line profile uses users provided Doppler HWHM, so you need to enter `Y`, `YES`, `Yes`, `y` or `yes` here and give a Doppler HWHM value.

Doppler line profile uses program calculated temperature-dependent Doppler HWHM. Program will not use the HWHM value from the input file.

###### Please type the correct Doppler HWHM choice 'Y' or 'N' into the input file.

When you choose Voigt profiles for calculating cross sections, you can provide a Doppler HWHM value (you need to enter `Y`, `YES`, `Yes`, `y` or `yes` here and give a Doppler HWHM value) or just let prgram to calculate temperature-dependent Doppler HWHM (you need to enter `N`, `NO`, `No`, `n`, `no` here).

###### Please type the correct Lorentzian HWHM choice 'Y' or 'N' into the input file.

When you choose Lorentzian or Voigt profiles for calculating cross sections, you can provide a Lorentzian HWHM value (you need to enter `Y`, `YES`, `Yes`, `y` or `yes` here and give a Lorentzian HWHM value) or just let prgram to calculate temperature and pressure dependent Lorentzian HWHM (you need to enter `N`, `NO`, `No`, `n`, `no` here).

###### Please add the name of the database 'ExoMol' or 'HITRAN' into the input file.

Make sure you have written the database name `ExoMol` or `HITRAN` in the input file.

###### The xxx boradening file does not exist.

Please check this xxx broadening file is in the correct folder, otherwise, you can use default broadening value to calculate cross sections.

###### The input file xxx does not exist.

Please check this xxx input file is in the correct folder and provide the correct input file path.

###### The file xxx is not a HITRAN2004 format data file.

If you use HITRAN database to do the calculations, please check your HITRAN data format. Only HITRAN2004 standard format can be read, otherwise, you need to convert your data format before using *PyExoCross*.

###### Please choose functions which you want to calculate.

Please write `1` after the function names which you want to calculate.

###### Please choose one from: 'Absoption' or 'Emission'.

Please enter `Absorption` or `Emission` in the input file.

###### Please choose wavenumber or wavelength and type in correct format: wn or wl.

Please enter `wn` (wavenumber) or `wl` (wavelength) in the input file.

###### Please choose line profile from the list.

Make sure you have written the line profile in the input file and please check your line profile name is correct.

###### Empty result with the input filter values. Please type new filter values in the input file.

Your uncertainty, threshold or quantum number filters may cause this problem so you have to reduce the filter range. You can increase uncertainty filter value, reduce threshold filter value or reduce quantum number filter range.

###### No uncertainties in states file. Please do not use uncertainty filter.

Please write `N`, `NO`, `No`, `n`, `no` after uncertainty filter.

###### concurrent.futures.process.BrokenProcessPool: A process in the process pool was terminated abruptly while the future was running or pending.

Each line list file is tooooo large when program running on many cores, so that the memory is not enough. Please reduce the `NCPUtrans`, `NCPUfiles` and/or `ChunkSize`. If you will use very large dataset (e.g. datasets of CaOH, NaOH, CH4, H2CS, VO), please set `NCPUfiles` as `1` and make `NCPUtrans` larger (try from 32).


###### ModuleNotFoundError: No module named 'distutils'.
Unfortunately, you're using the newly released Python, which removed distutils after it being deprecated since Python 3.10. Try to sideload distutils from a third party source (e.g. a system package, like apt install python3-distutils), or
downgrade to an older version of Python (3.10 or older). Recommend version 3.9 and 3.10.
