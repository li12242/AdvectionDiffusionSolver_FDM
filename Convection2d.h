//
//  ConvectionDrive.h
//  Project1
//
//  Created by li12242 on 16/10/3.
//  Copyright (c) 2016年 li12242. All rights reserved.
//

#ifndef Project1_ConvectionDrive_h
#define Project1_ConvectionDrive_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include "mkl.h"

typedef struct{
    int ndim2;  // num of points along x axis
    int ndim1;  // num of points along y axis
    double dx;  // distance between x points
    double dy;  // distance between y points
    double **x; // x coordinate
    double **y; // y coordinate
} structMesh;

typedef struct{
    double x0;
    double y0;
    double u;   // velocity
    double v;
    double Dx;
    double Dy;
} physics;

// ConvectionDrive2d.c

structMesh* MeshCreate();
void        MeshFree(structMesh *);

physics* PhysicsCreate();
void     PhysicsFree(physics *phys);

double** ConcentrationCreate(structMesh *mesh, physics *phys);
void     ConcentrationFree(double **);

void ExactSol(structMesh *mesh, physics *phys, double time, double **c_exa);
void NormErr (structMesh *mesh, physics *phys, double time, double **c,
              double *L2, double *Linf);


// ConvectionSolve2d.c
int ConvectionSolve2d(structMesh *mesh, physics *phys, double **c, double finalTime);

// ConvectionRHS2d.c
void ConvectionRHS2d_upwind (structMesh *mesh, physics *phys, double  **c, double **rhs);
void ConvectionRHS2d_central(structMesh *mesh, physics *phys, double  **c, double **rhs);

#endif //Project1_ConvectionDrive_h
