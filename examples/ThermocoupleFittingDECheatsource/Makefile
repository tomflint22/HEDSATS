# Makefile
# GNU makefile for example HEDSATS code solving 3D transient heat equation for double ellipsoidal heat source. Various Boundary conditions
# Questions/comments to thomas.flint@manchester.ac.uk (Tom Flint)

# includes
incdir = ../../HEDSATSlib

# compilers/flags
compiler = g++
flags = -O3 -I $(incdir) -std=c++11 -fopenmp -ftemplate-depth-1000 -Drestrict=


# dependencies
core = $(incdir)/double_exponential_transform_header.h \
       $(incdir)/Faddeeva.cpp \
       $(incdir)/Analytical_Functions.cpp

# the program
HScal: HeatSourceCalibrationWithDEC.cpp
	$(compiler) $(flags) $< -o $@


clean:
	rm -f HScal cooling.dat Output.dat
