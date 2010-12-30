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
OBJECTS      = calcrdg.o clearatm.o grillado.o imprime.o \
               islong.o petu.o poneres.o readdata.o setr.o  isclash.o
FILES        = calcrdg.cpp clearatm.cpp grillado.cpp \
               imprime.cpp isclash.cpp islong.cpp petu.cpp poneres.cpp readdata.cpp setr.cpp
HEADERS      = petu.h
LIBS         = -lprot-filer -lm -lgetopt_pp
COMBINATIONS_DEBUG_FLAGS = -DCOMBINATIONS_DEBUG
VERBOSE_FLAGS      = -DVERBOSE
DEBUG = -ggdb3
#LDFLAGS      = -L/usr/local/gromacs/lib   


CXXFLAGS     =-DMILI_NAMESPACE -Wall
#CPPFLAGS     = -I/usr/local/gromacs/include/gromacs	

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
gdb_petu:	petu


gdb_debug:	CXXFLAGS+= $(DEBUG)
gdb_debug:	debug


debug:  	CXXFLAGS+= $(COMBINATIONS_DEBUG_FLAGS)
debug:  	petu

verbose:  	CXXFLAGS+= $(VERBOSE_FLAGS)
verbose:  	CFLAGS+= $(VERBOSE_FLAGS)
verbose:  	petu


		
clean:
	    rm *.o  

