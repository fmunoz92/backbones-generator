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
OBJECTS      = calcrdg.o clearatm.o copymat.o imprime.o int2car.o islong.o petu.o \
               poneres.o readdata.o semilla.o setr.o writextc.o isclash.o
FILES        = calcrdg.c clearatm.c copymat.c imprime.c int2car.c isclash.c islong.c \
               petu.c poneres.c readdata.c semilla.c setr.c writextc.c 
HEADERS      = petu.h
LIBS         = -lxdrfile -lm -lpopt
#LDFLAGS      = -L/usr/local/gromacs/lib   

CFLAGS       = -O3
#CPPFLAGS     = -I/usr/local/gromacs/include/gromacs	
CC           = gcc
CXX          = g++
LD           = $(CXX)

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

petu:	        $(HEADERS) $(OBJECTS) 
		$(LD) $(LDFLAGS) -o $@ $(OBJECTS)   $(LIBS)


clean:
	    rm *.o  

