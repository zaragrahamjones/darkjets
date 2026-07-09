# PYTHIA configuration file.
# Generated on Wed Jul 8 19:16:36 JST 2026 with the following command:
# ./configure --with-fastjet3=/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local --with-rivet=/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local

# Install directory prefixes.
PREFIX_BIN=/Users/tvcr49@durham.ac.uk/code_software/pythia8315/bin
PREFIX_INCLUDE=/Users/tvcr49@durham.ac.uk/code_software/pythia8315/include
PREFIX_LIB=/Users/tvcr49@durham.ac.uk/code_software/pythia8315/lib
PREFIX_SHARE=/Users/tvcr49@durham.ac.uk/code_software/pythia8315/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_DTAGS=
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread
CXX_SHARED=-dynamiclib
CXX_SONAME=-Wl,-dylib_install_name,@rpath/
LIB_SUFFIX=.dylib
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local/bin/
FASTJET3_INCLUDE=-I/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local/include
FASTJET3_LIB=-L/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local/lib -Wl,-rpath,/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=true
RIVET_CONFIG=rivet-config
RIVET_BIN=/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local/bin/
RIVET_INCLUDE=-I/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local/include
RIVET_LIB=-L/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local/lib -Wl,-rpath,/Users/tvcr49@durham.ac.uk/code_software/rivetbootstrap/local/lib -lRivet

YODA_USE=false
YODA_CONFIG=
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

ROOT_USE=false
ROOT_CONFIG=
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
