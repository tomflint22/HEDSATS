# Makefile
# GNU makefile for example HEDSATS code solving 3D transient heat equation for double ellipsoidal heat source. Various Boundary conditions
# Questions/comments to thomas.flint@manchester.ac.uk (Tom Flint)

# includes
incdir = ../../HEDSATSlib

# compilers/flags
compiler = g++
flags = -O3 -I $(incdir)


# dependencies
core = $(incdir)/double_exponential_transform_header.h \
       $(incdir)/Faddeeva.cpp \
       $(incdir)/Analytical_Functions.cpp

# the program
DEC: DEC.cpp
	$(compiler) $(flags) $< -o $@


clean:
	rm -f DEC
