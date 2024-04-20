## Installation
### Dependencies
Cayley-RPMDrate depends on several other packages in order to provide its full 


functional capabilities. The following dependencies are required to compile and run RPMD:


Python: any version of Python3 is recommended


Numpy: version 1.5.0 or later is recommended


A standard Fortran 90/95 compiler


### Compiling from Source
Change the compiler and python choices in the Makefile and run `make` to compile RPMD.
### Generate PES file
run command `python generate.py` to generate NEW.so file prepared for RPMDrate. The pipnn file is specially optimised to be several times more efficient for RPMD calculations.
