PETSC_DIR = 
PETSC_ARCH = arch-linux-c-opt

FC = mpiifort
libs =

CFLAGS = ${PETSC_CC_INCLUDES}
FFLAGS = ${PETSC_FC_INCLUDES}

include ${PETSC_DIR}/lib/petsc/conf/variables

PETSC_INCLUDES = ${PETSC_FC_INCLUDES}


objects = FMM_VARS.o VARS.o DISTRIBUTOR.o READ_INP.o READ_INP_FMM.o REORDER_MESH.o\
	FORWARD_FMM.o ASSEMBLE_FMM.o BUILD_AMAP.o JACOBIAN_FMM.o   \
	JACOBIAN.o mat_inv.o FORWARD.o ASSEMBLE.o BuildMesh.o \
	DD_OPT_SLAVE.o SLAVE.o SLAVE_FMM.o\
	BUILD_MODEL_CONSTRAINTS.o OBJ.o\
	REPORT.o REPORT_FMM.o Output.o Output_FMM.o MASTER.o MASTER_FMM.o  \
	Multi_Forward.o INVERT.o JOINT_INVERT.o FMM_MAIN.o\
	DD_OPT.o V_ANALYTIC.o E4D.o

e4d : $(objects)
	$(FC) $(objects) -o e4d ${PETSC_LIB}

E4D.o : E4D.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 E4D.F90  -D resmode -o E4D.o ${PETSC_FC_INCLUDES}

VARS.o : VARS.F90
	$(FC)  -c -heap-arrays -Ofast -real-size 64 VARS.F90 -o VARS.o ${PETSC_FC_INCLUDES}

MASTER.o : MASTER.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 MASTER.F90 -o MASTER.o ${PETSC_FC_INCLUDES} 

SLAVE.o : SLAVE.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 SLAVE.F90 -D resmode -o SLAVE.o ${PETSC_FC_INCLUDES} 

READ_INP.o : READ_INP.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 READ_INP.F90 -D resmode -D sip -o READ_INP.o ${PETSC_INCLUDES}

BUILD_AMAP.o : BUILD_AMAP.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 BUILD_AMAP.F90 -o BUILD_AMAP.o ${PETSC_FC_INCLUDES}   

mat_inv.o : mat_inv.F
	$(FC) -c -heap-arrays -Ofast -real-size 64 mat_inv.F -o mat_inv.o ${PETSC_FC_INCLUDES}

IMI.o : IMI.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 IMI.F90 -o IMI.o ${PETSC_FC_INCLUDES}

FORWARD.o : FORWARD.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 FORWARD.F90 -o FORWARD.o ${PETSC_FC_INCLUDES}

Multi_Forward.o : Multi_Forward.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 Multi_Forward.F90 -o Multi_Forward.o ${PETSC_FC_INCLUDES}

ASSEMBLE.o : ASSEMBLE.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 ASSEMBLE.F90 -o ASSEMBLE.o ${PETSC_FC_INCLUDES}

Output.o : Output.f90
	$(FC) -c -heap-arrays -Ofast -real-size 64 Output.f90 -o Output.o ${PETSC_INCLUDES}

BUILD_MODEL_CONSTRAINTS.o : BUILD_MODEL_CONSTRAINTS.f90
	$(FC) -c -heap-arrays -Ofast -real-size 64 BUILD_MODEL_CONSTRAINTS.f90 -o BUILD_MODEL_CONSTRAINTS.o ${PETSC_INCLUDES}

OBJ.o : OBJ.f90
	$(FC) -c -heap-arrays -Ofast -real-size 64 OBJ.f90 -o OBJ.o ${PETSC_INCLUDES}

REPORT.o : REPORT.f90
	$(FC) -c -heap-arrays -Ofast -real-size 64 REPORT.f90 -o REPORT.o ${PETSC_INCLUDES}

JACOBIAN.o : JACOBIAN.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 JACOBIAN.F90 -o JACOBIAN.o ${PETSC_INCLUDES}

INVERT.o : INVERT.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 INVERT.F90 -o INVERT.o ${PETSC_INCLUDES}

JOINT_INVERT.o : JOINT_INVERT.F90
	 $(FC) -c -heap-arrays -Ofast -real-size 64 JOINT_INVERT.F90 -o JOINT_INVERT.o ${PETSC_INCLUDES}

BuildMesh.o : BuildMesh.f90
	$(FC) -c -heap-arrays -Ofast -real-size 64 BuildMesh.f90 -o BuildMesh.o ${PETSC_INCLUDES}

V_ANALYTIC.o : V_ANALYTIC.f90
	$(FC) -c -heap-arrays -Ofast -real-size 64 V_ANALYTIC.f90 -o V_ANALYTIC.o ${PETSC_INCLUDES}

DD_OPT.o : DD_OPT.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 DD_OPT.F90 -o DD_OPT.o ${PETSC_INCLUDES}

DD_OPT_SLAVE.o : DD_OPT_SLAVE.F90
	$(FC)  -c -heap-arrays -Ofast -real-size 64 DD_OPT_SLAVE.F90 -o DD_OPT_SLAVE.o ${PETSC_INCLUDES}

FMM_VARS.o : FMM_VARS.F90
	   $(FC)  -c -heap-arrays -Ofast -real-size 64 FMM_VARS.F90 -o FMM_VARS.o ${PETSC_INCLUDES}

DISTRIBUTOR.o : DISTRIBUTOR.F90
	   $(FC)  -c -heap-arrays -Ofast -real-size 64 DISTRIBUTOR.F90 -o DISTRIBUTOR.o ${PETSC_INCLUDES}

FMM_MAIN.o : FMM_MAIN.F90
	   $(FC)  -c -heap-arrays -Ofast -real-size 64 FMM_MAIN.F90 -o FMM_MAIN.o ${PETSC_INCLUDES}

FORWARD_FMM.o : FORWARD_FMM.F90
	   $(FC)  -c -heap-arrays -Ofast -real-size 64 FORWARD_FMM.F90 -o FORWARD_FMM.o ${PETSC_INCLUDES}

MASTER_FMM.o : MASTER_FMM.F90
	   $(FC)  -c -heap-arrays -Ofast -real-size 64 MASTER_FMM.F90 -o MASTER_FMM.o ${PETSC_INCLUDES}

SLAVE_FMM.o : SLAVE_FMM.F90
	   $(FC)  -c -heap-arrays -Ofast -real-size 64 SLAVE_FMM.F90 -o SLAVE_FMM.o ${PETSC_INCLUDES}

ASSEMBLE_FMM.o : ASSEMBLE_FMM.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 ASSEMBLE_FMM.F90 -o ASSEMBLE_FMM.o ${PETSC_FC_INCLUDES}

Output_FMM.o : Output_FMM.f90
	$(FC) -c -heap-arrays -Ofast -real-size 64 Output_FMM.f90 -o Output_FMM.o ${PETSC_INCLUDES}

REPORT_FMM.o : REPORT_FMM.f90
	   $(FC)  -c -heap-arrays -Ofast -real-size 64 REPORT_FMM.f90 -o REPORT_FMM.o ${PETSC_INCLUDES}

READ_INP_FMM.o : READ_INP_FMM.F90
	   $(FC)  -c -heap-arrays -Ofast -real-size 64 READ_INP_FMM.F90 -o READ_INP_FMM.o ${PETSC_INCLUDES}
	   
JACOBIAN_FMM.o : JACOBIAN_FMM.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 JACOBIAN_FMM.F90 -o JACOBIAN_FMM.o ${PETSC_INCLUDES}

MULTI_FORWARD.o : Multi_Forward.F90
	$(FC) -c -heap-arrays -Ofast -real-size 64 Multi_Forward.F90 -o MULTI_FORWARD.o ${PETSC_INCLUDES}

REORDER_MESH.o : REORDER_MESH.f90
	$(FC) -c -heap-arrays -Ofast -real-size 64 REORDER_MESH.f90 -o REORDER_MESH.o ${PETSC_INCLUDES}

.PHONY: clean
clean :
	rm -f *.o *.mod e4d *~
