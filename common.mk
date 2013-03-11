#
# Definitions common to all make files.
#

HDF5INCLUDEDIR ?= ../../seymour/dist/common/include
HDF5LIBDIR     ?= ../../seymour/dist/common/lib/

INCLUDEDIRS = -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR)

HDF5LIB    ?= hdf5
HDF5LIBCPP = hdf5_cpp

GCCOPTS = -O3 -Wno-div-by-zero

CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS)
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  
CPP = g++

ifeq ($(shell gcc -dumpversion | awk -F '.' '$$1*100+$$2>404{print "yes"}'),yes)
    CPPOPTS += -fpermissive
endif
ifneq ($(shell uname -s),Darwin)
    STATIC   = -static
endif
