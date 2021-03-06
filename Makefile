EXE			+= main
SOURCES	+= main.cpp mesh.cpp residue.cpp
OBJECTS	+= main.o conditions.o mesh.o residue.o
#PETSC_ARCH=arch-linux2-cxx-debug
#PETSC_ARCH =arch-linux2-cxx-opt-cluster
PETSC_ARCH = arch-linux2-cxx-opt
#FEPIC_DIR = /home/stevens/Documentos/USP/Proyecto/cluster/espeisan/libs/FEPiCpp
FEPIC_DIR = /home/stevens/Documentos/repos/fepicpp/FEPiCpp
PETSC_DIR = /home/stevens/Documentos/USP/Proyecto/cluster/espeisan/libs/petsc/petsc-3.4.4
#PETSC_DIR = /home/stevens/Documentos/USP/Proyecto/cluster/espeisan/libs/petsc/petsc-3.6.2

CFLAGS   = -g -Wall -Wextra
FFLAGS   =
CPPFLAGS = -g -I. -I${FEPIC_DIR} $(FEP_INCLUDE) -Wall -Wextra -fopenmp -m64 -msse2 -L$(FEP_LIBS_DIR) -lfepic $(FEP_LDFLAGS) -I ead -D EAD_DEBUG
FPPFLAGS =

.PHONY: all clean

# tests some variables
ifeq "" "$(wildcard ${FEPIC_DIR})"
$(error variable FEPIC_DIR was not defined or is an invalid directory)
endif 

ifeq "" "$(wildcard ${PETSC_DIR})"
$(error variable PETSC_DIR was not defined or is an invalid directory)
endif 

ifeq "" "$(wildcard ${PETSC_DIR}/${PETSC_ARCH})"
$(error variable PETSC_ARCH was not defined or is an invalid directory)
endif 

ifeq "" "${PETSC_ARCH}"
$(error PETSC_ARCH is was properly defined)
endif


#onde é procurado as bibliotecas
#SLINKER		+= -Wl,-rpath,/home/felipe/Slibs/smesh -L/home/felipe/Slibs/smesh

# "variables" must be included twice
# for petsc-3.4.4
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
# for petsc-3.6.2
#include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules

include ${FEPIC_DIR}/conf/variables

PETSC_KSP_LIB += -L$(FEP_LIBS_DIR) -lfepic $(FEP_LDFLAGS)

all: main

main: ${OBJECTS} chkopts
	-${CLINKER} -o $(EXE) $(OBJECTS) ${PETSC_KSP_LIB}

conditions.o: conditions.cpp
	g++ -c $(FEP_INCLUDE) conditions.cpp -o conditions.o

clean::
	${RM} *.o

# for petsc-3.4.4
include ${PETSC_DIR}/conf/test
# for petsc-3.6.2
#include ${PETSC_DIR}/lib/petsc/conf/test
