# pivpr

Exploration of the Pressure Relaxation phenomena with Particle Image Velocimetry, or `pivpr` for short. 

This repository is for code and resources for my masters thesis in Aerospace Engineering at [Old Dominion University](https://www.odu.edu/mae). The code revolves around synthesis and display of test data taken of an axial wake vortex with Particle Image Velocitmetry (PIV) in the ODU Low Speed Wind Tunnel (LSWT).

The data is discussed in the context of the pressure relaxation phenomena which is theorized to allow the longevity of axial wake vortex flow structures. "Raw" vector data is calculated from raw stereo image data and instrument calibration parameters by comercial software using particle displacement tracking methods, and is stored in tabular `v3d` files. These v3d files are considered the starting point for this code base, though some additional code has been writen to demonstrate uncertainty principles behind correlation techniques used to calculate 3d vector fields from stereo frame stradled image pairs.

![A graph](/texdocs/figs/example_vortex_figs/example_TscatterTKE.png)
