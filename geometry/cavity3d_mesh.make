EXEC = cavity3d_mesh

HOST = gcc-openmp
#HOST = gcc
#HOST = intel
#HOST = intel-openmp

FMMBIE_INSTALL_DIR = $(PREFIX)
FMM_INSTALL_DIR = $(PREFIX_FMM)
LFMMLINKLIB = -lfmm3d
LFMMLINKLIB += -Wl,-rpath,$(FMM_INSTALL_DIR)
LLINKLIB = -lfmm3dbie
LLINKLIB += -Wl,-rpath,$(FMMBIE_INSTALL_DIR)

ifneq ($(OS),Windows_NT)
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        ifeq ($(PREFIX_FMM),)
            FMM_INSTALL_DIR=/usr/local/lib
        endif
        ifeq ($(PREFIX),)
            FMMBIE_INSTALL_DIR=/usr/local/lib
        endif
    endif
    ifeq ($(UNAME_S),Linux)
        ifeq ($(PREFIX_FMM),)
            FMM_INSTALL_DIR=${HOME}/lib
        endif
        ifeq ($(PREFIX),)
            FMMBIE_INSTALL_DIR=${HOME}/lib
        endif
    endif
endif

ifeq ($(HOST),gcc)
    FC=gfortran
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy
endif

ifeq ($(HOST),intel)
    FC=ifort
    FFLAGS= -O3 -fPIC -march=native
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort
    FFLAGS= -O3 -fPIC -march=native -qopenmp
endif

FEND = -L${FMMBIE_INSTALL_DIR} $(LLINKLIB)

.PHONY: all clean

OBJECTS = cavity3d_geom.o \
  cavity3d_mesh.o

%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@ $(FEND)

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@ $(FEND)

all: $(OBJECTS)
	mkdir -p meshes
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(FEND)
	./$(EXEC)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)
