include ${PETSC_DIR}/lib/petsc/conf/variables

FC =
FFLAGS =
libs =
PETSC_DIR = 
PETSC_INCLUDES = ${PETSC_FC_INCLUDES}

objects = utility.o vtk.o \
	VARS.o DISTRIBUTOR.o READ_INP.o REORDER_MESH.o \
	BUILD_AMAP.o \
	JACOBIAN.o mat_inv.o FORWARD.o ASSEMBLE.o BuildMesh.f90.o \
	DD_OPT_SLAVE.o SLAVE.o \
	BUILD_MODEL_CONSTRAINTS.o OBJ.o\
	REPORT.o Output.o MASTER.o  \
	Multi_Forward.o INVERT.o \
	DD_OPT.o V_ANALYTIC.o E4D.o

e4d : $(objects)
	$(FC) $(objects) -o e4d ${PETSC_LIB}

E4D.o : E4D.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c E4D.F90  -D resmode -o E4D.o ${PETSC_FC_INCLUDES}

VARS.o : VARS.F90
	$(FC)  $(FFLAGS) -r8 -heap-arrays -c VARS.F90 -o VARS.o ${PETSC_FC_INCLUDES}

MASTER.o : MASTER.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c MASTER.F90 -o MASTER.o ${PETSC_FC_INCLUDES} 

SLAVE.o : SLAVE.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c SLAVE.F90 -D resmode -o SLAVE.o ${PETSC_FC_INCLUDES} 

READ_INP.o : READ_INP.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c READ_INP.F90 -D resmode -D sip -o READ_INP.o ${PETSC_INCLUDES}

BUILD_AMAP.o : BUILD_AMAP.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c BUILD_AMAP.F90 -o BUILD_AMAP.o ${PETSC_FC_INCLUDES}   

mat_inv.o : mat_inv.F
	$(FC) $(FFLAGS) -r8 -heap-arrays -c mat_inv.F -o mat_inv.o ${PETSC_FC_INCLUDES}

FORWARD.o : FORWARD.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c FORWARD.F90 -o FORWARD.o ${PETSC_FC_INCLUDES}

Multi_Forward.o : Multi_Forward.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c Multi_Forward.F90 -o Multi_Forward.o ${PETSC_FC_INCLUDES}

ASSEMBLE.o : ASSEMBLE.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c ASSEMBLE.F90 -o ASSEMBLE.o ${PETSC_FC_INCLUDES}

Output.o : Output.f90
	$(FC) $(FFLAGS) -r8 -c Output.f90 -o Output.o ${PETSC_INCLUDES}

BUILD_MODEL_CONSTRAINTS.o : BUILD_MODEL_CONSTRAINTS.f90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c BUILD_MODEL_CONSTRAINTS.f90 -o BUILD_MODEL_CONSTRAINTS.o ${PETSC_INCLUDES}

OBJ.o : OBJ.f90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c OBJ.f90 -o OBJ.o ${PETSC_INCLUDES}

REPORT.o : REPORT.f90
	$(FC) $(FFLAGS) -r8 -c REPORT.f90 -o REPORT.o ${PETSC_INCLUDES}

JACOBIAN.o : JACOBIAN.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c JACOBIAN.F90 -o JACOBIAN.o ${PETSC_INCLUDES}

INVERT.o : INVERT.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c INVERT.F90 -o INVERT.o ${PETSC_INCLUDES}

BuildMesh.f90.o : BuildMesh.f90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c BuildMesh.f90 -o BuildMesh.f90.o ${PETSC_INCLUDES}

V_ANALYTIC.o : V_ANALYTIC.f90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c V_ANALYTIC.f90 -o V_ANALYTIC.o ${PETSC_INCLUDES}

DD_OPT.o : DD_OPT.F90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c DD_OPT.F90 -o DD_OPT.o ${PETSC_INCLUDES}

DD_OPT_SLAVE.o : DD_OPT_SLAVE.F90
	$(FC)  $(FFLAGS) -r8 -heap-arrays -c DD_OPT_SLAVE.F90 -o DD_OPT_SLAVE.o ${PETSC_INCLUDES}

DISTRIBUTOR.o : DISTRIBUTOR.F90
	   $(FC)  $(FFLAGS) -r8 -heap-arrays -c DISTRIBUTOR.F90 -o DISTRIBUTOR.o ${PETSC_INCLUDES}

REORDER_MESH.o : REORDER_MESH.f90
	$(FC) $(FFLAGS) -r8 -heap-arrays -c REORDER_MESH.f90 -o REORDER_MESH.o ${PETSC_INCLUDES}

utility.o : utility.F90
	$(FC) $(FFLAGS) -heap-arrays -c utility.F90 -o utility.o

vtk.o : vtk.F90
	$(FC) $(FFLAGS) -heap-arrays -c vtk.F90 -o vtk.o

.PHONY: clean
clean :
	rm -f *.o *.mod e4d *~
