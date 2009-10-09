# Generated automatically from Makefile.in by configure.
#
# This is a Gromacs 3.0 template makefile for your own utility programs.
#
# Copy this file to whatever directory you are using for your own
# software and add more targets like the template one below.
#
# If you are using gmake it is relatively straightforward to add
# an include based on environment variables (like previous Gromacs versions)
# to select compiler flags and stuff automatically, but below it is static:
#

# Variables set by the configuration script:
OBJECTS      = calcrdg.o clearatm.o copymat.o grillado.o imprime.o int2car.o islong.o petu.o \
               poneres.o readdata.o semilla.o setr.o writextc.o isclash.o
FILES        = calcrdg.c clearatm.c copymat.c grillado.cpp imprime.c int2car.c isclash.c islong.c \
               petu.cpp poneres.cpp readdata.cpp semilla.cpp setr.c writextc.c 
HEADERS      = petu.h
LIBS         = -lxdrfile -lm 
COMBINATIONS_DEBUG_FLAGS = -DCOMBINATIONS_DEBUG
VERBOSE_FLAGS      = -DVERBOSE
DEBUG = -ggdb3
#LDFLAGS      = -L/usr/local/gromacs/lib   

CFLAGS       = -O3 -Wall
CXXFLAGS     =  -Wall
#CPPFLAGS     = -I/usr/local/gromacs/include/gromacs	
CC           = gcc
CXX          = g++
LD           = $(CXX)

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

petu:	        $(HEADERS) $(OBJECTS) 
		$(LD) $(LDFLAGS) -o $@ $(OBJECTS)   $(LIBS)

# Si se habilita el modo debug se deshabilitan 
# todos los chequeos por lo que no se descarta ninguna 
# cadena.
gdb_petu:	CXXFLAGS+= $(DEBUG)
gdb_petu:	CFLAGS+= $(DEBUG)
gdb_petu:	petu


gdb_debug:	CXXFLAGS+= $(DEBUG)
gdb_debug:	CFLAGS+= $(DEBUG)
gdb_debug:	debug


debug:  	CXXFLAGS+= $(COMBINATIONS_DEBUG_FLAGS)
debug:  	petu

verbose:  	CXXFLAGS+= $(VERBOSE_FLAGS)
verbose:  	CFLAGS+= $(VERBOSE_FLAGS)
verbose:  	petu


		
clean:
	    rm *.o  

