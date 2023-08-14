files=rpmdrate/_main.pyf rpmdrate/_math.f90 rpmdrate/_surface.f90 rpmdrate/_main.f90 rpmdrate/blas_lapack.f90 -m rpmdrate._main
compiler= --compiler=intelem --fcompiler=intelem
python=python3

.PHONY: all build install test clean

all: build

build:

	$(python) -m numpy.f2py -c rpmdrate/_surface.f90 -m rpmdrate._surface 
	$(python) -m numpy.f2py -c $(files)  $(compiler)    

test:

	$(python) rpmdrate/test.py

clean:
	rm -rf build/
	rm -rf rpmdrate/*.pyc
	rm -rf rpmdrate/*.pyo
	rm -rf rpmdrate/*.pyd
	rm -rf rpmdrate/*.so 
	rm -rf rpmdrate/__pycache__/*pyc
