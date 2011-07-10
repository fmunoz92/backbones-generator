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
OBJECTS      = src/utils.o src/grillado.o src/poneres.o src/readdata.o src/tree_generator.o
FILES        = src/utils.cpp src/grillado.cpp src/tree_generator.cpp src/petu.cpp src/poneres.cpp src/readdata.cpp 
HEADERS      = src/petu.h
LIBS         = -lprot-filer -lm -lgetopt_pp
COMBINATIONS_DEBUG_FLAGS = -DCOMBINATIONS_DEBUG
VERBOSE_FLAGS      = -DVERBOSE
DEBUG = -ggdb3
#LDFLAGS      = -L/usr/local/gromacs/lib   


CXXFLAGS     =-DMILI_NAMESPACE -Wall -Wextra -pedantic
#CPPFLAGS     = -I/usr/local/gromacs/include/gromacs	

CXX          = g++
LD           = $(CXX)

# The real make targets - note that most make programs support
# the shortcut $^ instead of listing all object files a second
# time, but we cannot count on it...

.PHONY: clean install tools

petu:	$(HEADERS) $(OBJECTS) src/petu.o
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) src/petu.o  $(LIBS)

tools: 
	$(CXX) $(CXXFLAGS) $(LIBS) tools/read_xtc.cpp -o read_xtc
	$(CXX) $(CXXFLAGS) $(LIBS) tools/read_compressed.cpp -o read_compressed

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
	    rm -f src/*.o ; rm -f petu; rm -f read_xtc; rm -f read_compressed 

# Find all sources with available tests
TEST_SRCS := $(patsubst ./%, %, $(shell find|egrep "_test\.cpp$$") )
TEST_BINS := $(patsubst %.cpp, %, $(TEST_SRCS))

test: $(HEADERS) $(OBJECTS) 
	@# Run all tests
	for TEST in $(TEST_BINS); do \
		make "$$TEST"_run; \
		done

%_test: %_test.cpp
	@# Compile the source and link with the real .o
	g++ -Isrc $(OBJECTS) $@.cpp -o $@ \
		$(CXXFLAGS) $(LDFLAGS) \
		-lgtest_main -lgmock -lpthread -lprot-filer

%_test_run: %_test
	@# Zomg magic magic - this will make the test and then run it
	./$<
