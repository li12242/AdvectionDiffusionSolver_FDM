OBJS=ConvectionGPUDrive2d.o \
	utils.o 

CC	= icc
LD	= icpc
NVCC = nvcc
LIBS = -lcuda -lcudart

CFLAGS = -O3 -I$(CUDADIR)/include
CUDADIR= /usr/local/cuda/
LDFLAGS= -L$(CUDADIR)/lib


# types of files we are going to construct rules for
.SUFFIXES: .c .f .cu

# rule for .c files
.c.o:
	$(CC) $(CFLAGS) -o $*.o -c $*.c
.cu.o: 
	$(NVCC) -o $*.o -c $*.cu $(CFLAGS) 

main:ConvectionGPU2d ConvectionGPUTexture2d
	rm -r $(OBJS) ConvectionGPUSolve2d.o ConvectionGPUTextureSolve2d.o

ConvectionGPU2d:$(OBJS) ConvectionGPUSolve2d.o
	$(LD)  $(LDFLAGS) -o ConvectionGPU2d $(OBJS) ConvectionGPUSolve2d.o $(LIBS)
	

ConvectionGPUTexture2d:$(OBJS) ConvectionGPUTextureSolve2d.o
	$(LD)  $(LDFLAGS) -o ConvectionGPUTexture2d $(OBJS) ConvectionGPUTextureSolve2d.o $(LIBS)

run:
	./ConvectionGPU2d
	./ConvectionGPUTexture2d

clean:
	rm ConvectionGPU2d
	rm ConvectionGPUTexture2d
	rm *.txt