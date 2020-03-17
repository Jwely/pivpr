# pivpr

[![DOI](https://zenodo.org/badge/16003/Jwely/pivpr.svg)](https://zenodo.org/badge/latestdoi/16003/Jwely/pivpr)

Exploration of the Pressure Relaxation phenomena with Particle Image Velocimetry, or `pivpr` for short. 

## Abstract
Vortices are a common phenomenon in fluid flows that arise as kinetic 
	energy dissipates into heat via viscous interaction. They arise naturally 
	at large scales in the form of dust devils, tornadoes, and as a 
	counter-rotating vortex pair in the wake of 
	aircraft. It is important to understand the conditions leading to their 
	formation, their duration, and their dissipation in order to forecast or 
	prevent undesirable effects. Among these deleterious effects is a decrease 
	in safety of aircraft 
	operations in the wake of other aircraft, an extremely common situation at 
	airports around the world. A large number of mathematical models and 
	experimental data sets exists to help explain various aspects of axial wake 
	vortex behavior, but current models fail to explain why many vortices 
	remain tightly wound with slowly decaying azimuthal 
	velocities about their cores the length of time for which they 
	have been observed. The current study builds upon the theoretical work of 
	Ash, Zardadkhan and Zuckerwar, and tests specific attributes 
	of a turbulent axial vortex for agreement with non-equilibrium pressure 
	relaxation theory. This theory provides an exact solution to a modified 
	version of the Navier-Stokes equations for an 
	axial vortex, with a resulting velocity model that agrees with leading 
	empirical models. In the present investigation, axial wake vortices were 
	created with a bi-wing vortex 
	generator in a low speed wind tunnel, at free stream velocities between 15 
	and 33 $m/s$. Stereo particle image velocimetry was employed to 
	map three dimensional velocity vectors at positions between 5.4 and 10 
	chord lengths downstream of the vortex generator, and at a sampling rate of 
	1Hz for 200 seconds. A Reynolds time averaging approach was employed to 
	express instantaneous velocity measurements as localized mean and 
	fluctuating components and to study turbulent structures within the 
	vortices. Periodicity in turbulent energy and Reynolds stress structures 
	was observed by comparing vortex velocity fields normalized by age, based 
	on free stream velocity and downstream distance. The cores of these 
	vortices appeared to periodically ingest turbulent energy and compress it 
	into approximately one half of local core radii. The cyclical ingestion of 
	turbulence was shown to have the effect of tightening the core radius in 
	the wake of the vortex generator center body. If this phenomenon persists 
	for the life of the vortex, it could provide an explanation for the 
	longevity of the azimuthal velocity component, as observed in natural wake 
	vortices.
	
## About
This repository is for code and resources for my masters thesis in Aerospace Engineering at [Old Dominion University](https://www.odu.edu/mae). The code revolves around synthesis and display of test data taken of an axial wake vortex with Particle Image Velocitmetry (PIV) in the ODU Low Speed Wind Tunnel (LSWT).

The data is discussed in the context of the pressure relaxation phenomena which is theorized to allow the longevity of axial wake vortex flow structures. "Raw" vector data is calculated from raw stereo image data and instrument calibration parameters by comercial software using particle displacement tracking methods, and is stored in tabular `v3d` files. These v3d files are considered the starting point for this code base, though some additional code has been writen to demonstrate uncertainty principles behind correlation techniques used to calculate 3d vector fields from stereo frame stradled image pairs.

## Dataset and document access
The full 3d vector dataset and working drafts of the thesis can be downloaded from the [release page](https://github.com/Jwely/pivpr/releases). The vector dataset is 800MB compressed, but unzips to 12GB. The raw images used to synthesize the vector fields are not available for download due to the large file size, but may be supplied upon request.

## Requirements
Statistics sets and plots can be generated from the vector set using controller scripts which interface with a small axial vortex API. Runing this code requires:
 * [Python 2.7 +](https://www.python.org/downloads/)
 * [Microsoft visual C++ compiler (windows)](https://www.microsoft.com/en-us/download/details.aspx?id=44266)

And, with easy installation using pip:
 * [`numpy` with the math kernel library](http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy)
 * [`pandas`](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pandas)
 * [`matplotlib`](http://www.lfd.uci.edu/~gohlke/pythonlibs/#matplotlib)
 * [`libtiff`](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pylibtiff)
 * [`scipy`](http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy)

## Directory Structure
#### py
Living version of the codebase. Performs analysis and dumps files to texdocs. Built to make the analysis as repeatable as posisble, with full exposure to the methods used as part of this study.
#### matlab
Obsolete code base started in matlab, kept for good record keeping only.
#### texdocs
TeX and resources used to produce PDF formatted publication document. Several subdirectories are automatically filled by `py/controler` scripts with TeX files, figure jpgs, and tables converted from .csv format.
#### data_full
Folder to place full set of processed 3d vector files. Controllers in the `py\controler` package will work if data is correctly placed in this folder.
#### data_test
Small limited subset of vector files initially used for testing and development.

## Notes
The intention was to eventually use [OpenPIV](https://github.com/OpenPIV/openpiv-python) to perform computation, but that repository at the time of this research still needed a lot of work. This approach was abandoned, but I encourage future users of PIV and python to check it out.
