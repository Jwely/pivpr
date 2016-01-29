# pivpr

Exploration of the Pressure Relaxation phenomena with Particle Image Velocimetry, or `pivpr` for short. 

This repository is for code and resources for my masters thesis in Aerospace Engineering at [Old Dominion University](https://www.odu.edu/mae). The code revolves around synthesis and display of test data taken of an axial wake vortex with Particle Image Velocitmetry (PIV) in the ODU Low Speed Wind Tunnel (LSWT).

The data is discussed in the context of the pressure relaxation phenomena which is theorized to allow the longevity of axial wake vortex flow structures. "Raw" vector data is calculated from raw stereo image data and instrument calibration parameters by comercial software using particle displacement tracking methods, and is stored in tabular `v3d` files. These v3d files are considered the starting point for this code base, though some additional code has been writen to demonstrate uncertainty principles behind correlation techniques used to calculate 3d vector fields from stereo frame stradled image pairs.

## Dataset access
- The full 3d vector dataset for all experiments is about 13GB once decompressed and can be downloaded from the [release page](https://github.com/Jwely/pivpr/releases). The raw images used to synthesize the vector fields are not available for download due to the large file size.

- A `.7z` file with all the graphs anyone could possibly want is also available on the releases page.

## Directory Structure
#### py
Living version of the codebase. performs analysis and dumps files to texdocs. Many sub packages.
#### matlab
Outdated codebase, kept for good record keeping.
#### texdocs
Tex and tex resources used to produce PDF and HTML formatted research document.
#### data_full
Folder to place full set of processed 3d vector files. Controllers in the `py` package will work with data placed in this folder.
#### data_test
Small limited subset of vector files initially used for testing and development.

#### A graph! yay graphs!
![A graph](/texdocs/figs/example_vortex_figs/example_TscatterTKE.png)
