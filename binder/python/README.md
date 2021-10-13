# Custom Python Libraries

This directory houses libraries that we have written. (? I assume Anthony wrote dataPython but I should check)

## Libraries

`dataPython.py` is for importing data and uncertainties from text files to numpy arrays.

`fitting_NGC5533.py` expands upon `NGC5533_functions.py` with functions frequently used in our fitting routines.

`NGC5533_functions.py` are general equations for the shapes of components of the rotation curve for the galaxy NGC5533. It also allows saving and loading of .hdf5 files, and handles this in a way designed to work with the other functions.

`noordermeer.py` loads data traced from graphs in [Noordermeer 2018](https://doi.org/10.1111/j.1365-2966.2008.12837.x).
