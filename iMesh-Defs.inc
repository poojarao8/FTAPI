IMESH_CXX= /nfs/scratch/jmelvin/lib/openmpi/bin/mpicxx
IMESH_CXXFLAGS = $(CXXFLAGS) -g
IMESH_CFLAGS = $(CFLAGS) -g
IMESH_LDFLAGS = -L$(FT_DIR)/lib/x86_64   -L/nfs/scratch/lib/openmpi/include/lib 
IMESH_LIBS = -lFronTier

IMESH_F77 = /nfs/scratch/jmelvin/lib/openmpi/bin/mpif77
IMESH_F77_LIBS= 

IMESH_INCLUDES = -I$(FT_DIR)/include/FronTier  -I$(FT_DIR)/include/FronTier/intfc

IMESH_LIBS = -lFronTier -lm  -lmpi -lmpi_cxx  

IMESH_INCLUDEDIR = include/FronTier/intfc include/FronTier/front include/FronTier/util

FT_DIR=/nfs/scratch/jmelvin/software/FTAPI
