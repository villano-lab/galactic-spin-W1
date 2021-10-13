This is not thorough documentation on the functions, but I figured a few notes might be handy to have around

### Dependencies
The following dependencies are required:
- numpy
- scipy
The following dependency is optional, but highly recommended (several features will not work without it):
- h5py

### Function Naming Convention
All four contributions to the velocity have a prefix. 
These are 'h' (halo), 'bh' (black hole), 'd' (disk), and 'b' (bulge).
Any functions from one of these components will have the prefix. (Ex: `b_gammafunc`) </br>
All velocity functions are abbreviated to 'v'. (Ex: `bh_v`)
The total velocity function is simply titled `v`; it has no prefix.

### Function Usage
With the exception of `savedata` and `loaddata`, all functions have defaults such that only one input is required, but more may be allowed.
For instance, `bh_v(r)` will assume the measured value for the black hole mass, and `bh_v(r,1)` will assign that mass a value of 1.

All velocities functions have two key-word arguments that default to false: `save` and `load`. </br>
If you set `save=True`, the calculations will be saved to a file so that they can be loaded later.
If there is already existing data for your set variables, it will do all of the calculations and add those to whatever existing data there is while removing redundancies (points with the same inputs). </br>
If you set `load=True`, existing calculations for your parameters will be loaded.
These calculations are then splined and your input values will be plugged into the spline to get your output values.
If no data can be loaded for your values, one of two things will happen:
- If the function's only parameter can be treated as a prefactor, it will scale the data for which that prefactor is one.
    If the data for the prefactor of one cannot be loaded, it will throw an error.
- If the function has parameters which cannot be treated as prefactors, it will run the calculation fully and save the data for future use.
The functions are not written to handle both saving and loading.
They won't throw an error at setting both true, but only one will happen (which one is listed first depends on the function due to some issues with structuring error handling).
However, since load will redirect to saving if necessary, it is best to choose whichever you want to prioritize.</br>
All velocity functions have defaults set that lead them to use the expected hdf5 structure for saving and loading.
If you want to use a different structure, kwargs for the `savedata` and `loaddata` can be passed through the velocity functions.

### Expected hdf5 Structure
The `savedata` and `loaddata` functions use h5py to read and write to .hdf5 files.
Functions which save or load data have defaults set to use the expected structure for ease of use.
This is as follows:
- The file used for storing calculations is in the current directory. You can change this with the `path` variable.
- The name of this file is based on the component ('total.hdf5','disk.hdf5','bulge.hdf5','blackhole.hdf5').
- Currently, each file has only one group with the same name (without the hdf5). This structure is due to file size limitations.
- Within that file and group, each unique set of parameters will be a single dataset.
- Each dataset will contain two numpy arrays, the zeroth corresponding to x-data and the first corresponding to calculated y-values.