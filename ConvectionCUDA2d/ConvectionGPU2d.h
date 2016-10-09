//
//  ConvectionDrive.h
//  Project1
//
//  Created by li12242 on 16/10/3.
//  Copyright (c) 2016年 li12242. All rights reserved.
//

#ifndef ConvectionCUDADrive_h
#define ConvectionCUDADrive_h

#define NAMELEN 256
#define np 401 // number of points along single axis

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
// #inclu de "mkl.h"

typedef struct{
    int ndim2;  // num of points on 2nd dimension (y)
    int ndim1;  // num of points on 1st dimension (x)
    float dx;   // distance between x points
    float dy;   // distance between y points
    float **x;  // x coordinate matrix
    float **y;  // y coordinate matrix
} structMesh;

typedef struct{
    float x0;
    float y0;
    float u;   // velocity
    float v;   // velocity
    float Dx;
    float Dy;
} physics;

// ConvectionDrive2d.c

structMesh* MeshCreate();
void        MeshFree(structMesh *);

physics* PhysicsCreate();
void     PhysicsFree(physics *phys);

float** ConcentrationCreate(structMesh *mesh, physics *phys);
void     ConcentrationFree(float **);

void ExactSol(structMesh *mesh, physics *phys, float time, float **c_exa);
void NormErr (structMesh *mesh, physics *phys, float time, float **c,
              float *L2, float *Linf);


// ConvectionSolve2d.c
void ConvectionGPUSolve2d(structMesh *mesh, physics *phys, 
    float **c, float finalTime);

#endif //ConvectionCUDADrive_h
