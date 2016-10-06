Src=ConvectionParDrive2d.c \
	ConvectionParallel2d.h \
	ConvectionRHS2d.c \
	ConvectionSolve2d.c \
	../Utilities/utils.c \
	../Utilities/utils.h \

ConvectionParallel2d:
	mpicc ${Src} -o $@

# N is the num of processes
# Usages: 
#	make run N=3
run:
	mpirun -n ${N} -host localhost ./ConvectionParallel2d

clean:
	rm *.txt
	rm ConvectionParallel2d
	
