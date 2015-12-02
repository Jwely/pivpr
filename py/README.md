# thesis-pivpr/py
Code for synthesis of PIV data for demonstration of pressure relaxation phenomena. This is the second and better version of the codebase rewritten in python.

Thesis data was taken of an axial wake vortex in the ODU low speed wind tunnel with a stereo particle image velocimetry system. 
This code is built to synthesize this raw cartesian PIV imagery into a 3d velocity vector field in cylindrical coordinates about
the center of the vortex.

## to-do

- [ ] create controler scripts to construct AxialVortex instances and generate output plots procedurally (image format?)
- [ ] set plotting methods to auto-bound image extents as function of distance from vortex core
- [ ] generalize plotting methods to return figure handles so a parent function can nest subfigures easily
- [ ] move vortex specific methods from MeanVecFieldCylindrical to AxialVortex
- [ ] finish scatter plotting methods
