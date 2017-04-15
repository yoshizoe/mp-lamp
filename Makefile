## Set USE_LP to 1 to enable linear programming stuff.
USE_LP=0

## Set LINK_RELEASE_STATICALLY to 0 or 1 (default) to disable/enable
## static linking of the executable in release mode.
## On OS X, this is unsupported and will be silently disabled.
LINK_RELEASE_STATICALLY=1

## Set STATE_VAR_BYTES to 1, 2 or 4 to control how large a state
## variable is. There should be no need to tweak this manually.
STATE_VAR_BYTES=1

## On a supported operating system, there should be no need to override
## the OS setting. If the provided code does not work even though your
## operating system is a supported one, please report this as a bug.
OS=auto


HEADERS = \
	src/database.h \
	src/dfs.h \
	src/graph.h \
	src/hash.h \
	src/lamp_graph.h \
	src/lamp.h \
	src/sorted_itemset.h \
	src/table_vba.h \
	src/table.h \
	src/timer.h \
	src/topk.h \
	src/utils.h \
	src/variable_bitset_array.h \
	src/variable_length_array.h

HEADERS += \
	mp-src/DTD.h \
	mp-src/Log.h \
	mp-src/SignificantSetResults.h \
	mp-src/StealState.h \
	mp-src/mp_dfs.h


SHELL = /bin/bash

SOURCES = mp-src/mp-lamp.cc $(HEADERS:%.h=%.cc)
TARGET  = mp-lamp



OBJECT_SUFFIX_RELEASE = .$(STATE_VAR_BYTES)
TARGET_SUFFIX_RELEASE = -$(STATE_VAR_BYTES)
OBJECT_SUFFIX_DEBUG   = .$(STATE_VAR_BYTES).debug
TARGET_SUFFIX_DEBUG   = -$(STATE_VAR_BYTES)-debug
OBJECT_SUFFIX_PROFILE = .$(STATE_VAR_BYTES).profile
TARGET_SUFFIX_PROFILE = -$(STATE_VAR_BYTES)-profile

OBJECT_SUFFIX_MPI_PROFILE = .$(STATE_VAR_BYTES).mpiprofile
TARGET_SUFFIX_MPI_PROFILE = -$(STATE_VAR_BYTES)-mpiprofile

OBJECTS_RELEASE = $(SOURCES:%.cc=.obj/%$(OBJECT_SUFFIX_RELEASE).o)
TARGET_RELEASE  = $(TARGET)$(TARGET_SUFFIX_RELEASE)

OBJECTS_DEBUG   = $(SOURCES:%.cc=.obj/%$(OBJECT_SUFFIX_DEBUG).o)
TARGET_DEBUG    = $(TARGET)$(TARGET_SUFFIX_DEBUG)

OBJECTS_PROFILE = $(SOURCES:%.cc=.obj/%$(OBJECT_SUFFIX_PROFILE).o)
TARGET_PROFILE  = $(TARGET)$(TARGET_SUFFIX_PROFILE)

OBJECTS_MPI_PROFILE = $(SOURCES:%.cc=.obj/%$(OBJECT_SUFFIX_MPI_PROFILE).o)
TARGET_MPI_PROFILE  = $(TARGET)$(TARGET_SUFFIX_MPI_PROFILE)

CC     = mpic++ #g++
DEPEND = mpic++ -MM 

## CCOPT, LINKOPT, POSTLINKOPT are options for compiler and linker
## that are used for all three targets (release, debug, and profile).
## (POSTLINKOPT are options that appear *after* all object files.)

CCOPT = -Iext
CCOPT += -g
CCOPT += -std=c++11
CCOPT += #-m32
CCOPT += -Wall -W -Wno-unused-parameter -Wno-variadic-macros -Wno-sign-compare -Wno-deprecated -pedantic -Werror -DSTATE_VAR_BYTES=$(STATE_VAR_BYTES)
CCOPT += -L/usr/lib
CCOPT += -I${JEMALLOC_PATH}/include -L${JEMALLOC_PATH}/lib -Wl,-rpath,${JEMALLOC_PATH}/lib

## The following lines contain workarounds for bugs when
## cross-compiling to 64 bit on 32-bit systems using gcc 4.4 or gcc
## 4.5 in some Ubuntu releases. (We don't usually cross-compile to
## 64-bit, but in some cases we do; e.g. we did for the IPC.) See
## http://stackoverflow.com/questions/4643197/missing-include-bits-cconfig-h-when-cross-compiling-64-bit-program-on-32-bit.

#HAVE_GCC_4_4 := $(shell expr "$$(gcc -dumpversion)" : \\\(4\.4\.\\\))
#HAVE_GCC_4_5 := $(shell expr "$$(gcc -dumpversion)" : \\\(4\.5\.\\\))

#ifdef HAVE_GCC_4_4
#    CCOPT += -I/usr/include/c++/4.4/i686-linux-gnu
#endif

#ifdef HAVE_GCC_4_5
#    CCOPT += -I/usr/include/c++/4.5/i686-linux-gnu
#endif

LINKOPT  = -g
LINKOPT += -L/usr/include/x86_64-linux-gnu/c++/4.8/bits/ #-m32 
LINKOPT += -L/usr/lib #-lmpi_cxx  # -lmpi -lhwloc -Wl,-rpath

POSTLINKOPT =

## Additional specialized options for the various targets follow.
## In release mode, we link statically since this makes it more likely
## that local compiles will work on the various grids (gkigrid, Black
## Forest Grid, maia).
##
## NOTE: This precludes some uses of exceptions.
##        For details, see man gcc on -static-libgcc.

CCOPT_RELEASE  = -O3 -DNDEBUG -fomit-frame-pointer
CCOPT_DEBUG    = -O3 -DDEBUG
CCOPT_PROFILE  = -O3 -pg 
CCOPT_MPI_PROFILE = -O3 -L/usr/local/lib -lmpiP -lm -lunwind

LINKOPT_RELEASE =
LINKOPT_DEBUG    =
LINKOPT_PROFILE  = -pg
LINKOPT_MPI_PROFILE = -L/usr/local/lib -lmpiP -lm -lunwind

#POSTLINKOPT_RELEASE = -ljemalloc
#POSTLINKOPT_DEBUG   = -ljemalloc
#POSTLINKOPT_PROFILE = -ljemalloc

ifeq ($(LINK_RELEASE_STATICALLY), 1)
LINKOPT_RELEASE += #-static -static-libgcc
endif

ifeq ($(OS), linux)
ifeq ($(LINK_RELEASE_STATICALLY), 0)
POSTLINKOPT_RELEASE += -lrt
else
POSTLINKOPT_RELEASE += #-Wl,-Bstatic -lrt
endif
POSTLINKOPT_DEBUG  += -lrt
POSTLINKOPT_PROFILE += -lrt
# POSTLINKOPT_PROFILE += -lc_p
endif

POSTLINKOPT_MPI_PROFILE = -lmpiP -lm -lunwind -lrt -ljemalloc


## Define the default target up here so that the LP stuff below
## doesn't define a default target.

default: release



all: release debug profile 

## Build rules for the release target follow.

release: $(TARGET_RELEASE)

$(TARGET_RELEASE): $(OBJECTS_RELEASE)
	$(CC) $(LINKOPT) $(LINKOPT_RELEASE) $(OBJECTS_RELEASE) $(POSTLINKOPT) $(POSTLINKOPT_RELEASE) -o $(TARGET_RELEASE)

$(OBJECTS_RELEASE): .obj/%$(OBJECT_SUFFIX_RELEASE).o: %.cc
	@mkdir -p $$(dirname $@)
	$(CC) $(CCOPT) $(CCOPT_RELEASE) -c $< -o $@

## Build rules for the debug target follow.

debug: $(TARGET_DEBUG)

$(TARGET_DEBUG): $(OBJECTS_DEBUG)
	$(CC) $(LINKOPT) $(LINKOPT_DEBUG) $(OBJECTS_DEBUG) $(POSTLINKOPT) $(POSTLINKOPT_DEBUG) -o $(TARGET_DEBUG)

$(OBJECTS_DEBUG): .obj/%$(OBJECT_SUFFIX_DEBUG).o: %.cc
	@mkdir -p $$(dirname $@)
	$(CC) $(CCOPT) $(CCOPT_DEBUG) -c $< -o $@

## Build rules for the profile target follow.

profile: $(TARGET_PROFILE)

$(TARGET_PROFILE): $(OBJECTS_PROFILE)
	$(CC) $(LINKOPT) $(LINKOPT_PROFILE) $(OBJECTS_PROFILE) $(POSTLINKOPT) $(POSTLINKOPT_PROFILE) -o $(TARGET_PROFILE)

$(OBJECTS_PROFILE): .obj/%$(OBJECT_SUFFIX_PROFILE).o: %.cc
	@mkdir -p $$(dirname $@)
	$(CC) $(CCOPT) $(CCOPT_PROFILE) -c $< -o $@


## Build rules for the MPIPROFILE target follow.

mpi_profile: $(TARGET_MPI_PROFILE)

$(TARGET_MPI_PROFILE): $(OBJECTS_MPI_PROFILE)
	$(CC) $(LINKOPT) $(LINKOPT_MPI_PROFILE) $(OBJECTS_MPI_PROFILE) $(POSTLINKOPT) $(POSTLINKOPT_MPI_PROFILE) -o $(TARGET_MPI_PROFILE)

$(OBJECTS_MPI_PROFILE): .obj/%$(OBJECT_SUFFIX_MPI_PROFILE).o: %.cc
	@mkdir -p $$(dirname $@)
	$(CC) $(CCOPT) $(CCOPT_MPI_PROFILE) -c $< -o $@

## Additional targets follow.

PROFILE: $(TARGET_PROFILE)
	./$(TARGET_PROFILE) $(ARGS_PROFILE)
	gprof $(TARGET_PROFILE) | (cleanup-profile 2> /dev/null || cat) > PROFILE

clean:
	rm -rf .obj
	rm -f *~ *.pyc
	rm -f Makefile.depend gmon.out PROFILE core
	rm -f sas_plan

distclean: clean
	rm -f $(TARGET_RELEASE) $(TARGET_DEBUG) $(TARGET_PROFILE)

## NOTE: If we just call gcc -MM on a source file that lives within a
## subdirectory, it will strip the directory part in the output. Hence
## the for loop with the sed call.

Makefile.depend: $(SOURCES) $(HEADERS)
	rm -f Makefile.temp
	for source in $(SOURCES) ; do \
	    $(DEPEND) $$source > Makefile.temp0; \
	    objfile=$${source%%.cc}.o; \
	    sed -i -e "s@^[^:]*:@$$objfile:@" Makefile.temp0; \
	    cat Makefile.temp0 >> Makefile.temp; \
	done
	rm -f Makefile.temp0 Makefile.depend
	sed -e "s@\(.*\)\.o:\(.*\)@.obj/\1$(OBJECT_SUFFIX_RELEASE).o:\2@" Makefile.temp >> Makefile.temp0
	sed -e "s@\(.*\)\.o:\(.*\)@.obj/\1$(OBJECT_SUFFIX_DEBUG).o:\2@" Makefile.temp >> Makefile.temp0
	sed -e "s@\(.*\)\.o:\(.*\)@.obj/\1$(OBJECT_SUFFIX_PROFILE).o:\2@" Makefile.temp >> Makefile.temp0
	rm -f Makefile.temp
	sed -e "s@\(.*\)\.$(STATE_VAR_BYTES)\.\(.*\)@\1.1.\2@" Makefile.temp0 >> Makefile.depend
	sed -e "s@\(.*\)\.$(STATE_VAR_BYTES)\.\(.*\)@\1.2.\2@" Makefile.temp0 >> Makefile.depend
	sed -e "s@\(.*\)\.$(STATE_VAR_BYTES)\.\(.*\)@\1.4.\2@" Makefile.temp0 >> Makefile.depend
	rm -f Makefile.temp0

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
-include Makefile.depend
endif
endif

.PHONY: default all release debug profile clean distclean
