Code for synthesis of PIV data for demonstration of pressure relaxation phenomena. This is the second and better version of the codebase rewritten in python.

Thesis data was taken of an axial wake vortex in the ODU low speed wind tunnel with a stereo particle image velocimetry system. 
This code is built to synthesize this raw cartesian PIV imagery into a 3d velocity vector field in cylindrical coordinates about
the center of the vortex.

### controler
Contains controller scripts for the other python modules within this directory. These controllers do everything
from analyze data, create TeX figures and tables, and perform uncertainty analysis. 

### piv
Directly manipulates and manages PIV data and provides an API for creating plots and statistics. The structure is heirarchal, and looks like

```
VecFieldCartesian             # an individual vector field instance created from csv like `.v3d` files
MeanVecFieldCartesian         # averages of many VecFieldCartesian instances, gets flow properties in cartesian
    ↳ AxialVortex             # extends parent with cylindrical coordinate analysis methods and plotting methods
Experiment                    # class to manage one AxialVortex and anciliary experimental data as one unit
construct_axial_vortex        # constructs an AxialVortex instance from directory of data, or loads previous
construct_experiment          # pairs an AxialVortex with anciliary data parsed from csv files

shorthand_to_tex              # Translates shorthand notations used to code to TeX formatted math symbols

\pickles                      # holds binaries of AxialVortex instances for quick lookup without recomputation
\dat                          # directory with anciliary data read into Experiments
```

### tex
Contains code specifically for creating TeX formatted figures and tables for use in the actual thesis document
generated with TeX. 

### uncertainty
Limited amount of code for synthesizing artificial PIV data used for monte carlo uncertainty analysis.

### utils
General utilities folder of short little useful functions. Includes a dependencies installer, a csv_to_tex converter, 
and a few other things.
