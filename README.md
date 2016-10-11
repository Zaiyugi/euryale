## Final Project: SPH2D with material editing and FLIP-style color grids
### Author: Zachary Shore
### Class: CPSC-881: Real-time Fluid Simulation

### Description

*  OpenImageIO: [Download](https://github.com/OpenImageIO/oiio)  

### Building
*  make -- make all required libraries
*  make final\_sph2d -- makes the final\_sph2d executable; requires libraries to be built

### Executables
*  final\_sph2d: Distribute particles based on specified parameters  

### Directories
*  _colormaps/_  
   Contains colormap CSV files; these are used by the ColorGradient class
*  _include/_  
   Contains class header files;
*  _lib/_  
   Contains static library files; three are made during compilation: libOGL.a, libVR.a, libMath.a
*  _ogl/_  
   Contains OpenGL headers and source files.
*  _src/_  
   Contains some source files

### Namespaces
*  _lux_: Most of the more mathematical classes are under this; Vector, Tensor, Matrix, etc.
*  _ogl_: Everything in opengl/ uses this
*  _util_: Contains utility functions; Check Utility.h for definitions
*  _sim_: Anything directly related to the sim code: SPHEngine, SPHForce, Particle, etc.

Created: 2016-04-01 
Edited: 2016-04-14
 
