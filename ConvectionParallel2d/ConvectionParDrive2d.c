//
//  main.c
//  ConvectionParallel2d
//
//  Created by li12242 on 16/10/5.
//  Copyright (c) 2016年 li12242. All rights reserved.
//

#include "ConvectionParallel2d.h"
#define np 201

int main(int argc, char **argv){
    
    MPI_Init(&argc, &argv);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    structMesh *mesh = ParallelMeshCreate();
    physics    *phys = PhysicsCreate();
    double     **c   = ConcentrationCreate(mesh, phys);
    
    double finalTime = 2.5;
    
    PrintTotalError(mesh, phys, c, 0, procid);
    ConvectionSolve2d(mesh, phys, c, finalTime);
    PrintTotalError(mesh, phys, c, finalTime, procid);
    
    char filename[200];
    snprintf(filename, 200, "result-%d-%d.txt", procid, nprocs);
//    printf("procid=%d, filename=%s\n", procid, filename);
    SaveMatrix(filename, c, mesh->ndim1, mesh->ndim2);
    snprintf(filename, 200, "y-%d-%d.txt", procid, nprocs);
    SaveMatrix(filename, mesh->y, mesh->ndim1, mesh->ndim2);
    snprintf(filename, 200, "x-%d-%d.txt", procid, nprocs);
    SaveMatrix(filename, mesh->x, mesh->ndim1, mesh->ndim2);
    
    ConcentrationFree(c);
    ParallelMeshFree(mesh);
    PhysicsFree(phys);
    
    MPI_Finalize();
    return 0;
}

/* Print the Norm Error on screen */
void PrintTotalError(structMesh *mesh, physics *phys,
                     double **c, double time, int procid){
    double L2, Linf, gL2, gLinf;
    
    NormErr(mesh, phys, time, c, &L2, &Linf);
    MPI_Allreduce(&L2,   &gL2,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&Linf, &gLinf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //printf("procid=%d, Local Norm Error, L2=%f, Linf=%f\n",procid,L2,Linf);
    if (!procid) {
        printf("Total Norm Error, L2=%f, Linf=%f\n", gL2,gLinf);
    }
    return;
}

/**
 * @brief
 * Allocate Mesh strucutre for Parallel Mesh strucutre.
 *
 * @details
 * Divide rows of the mesh into each processes.
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh  | structMesh* | pointer to the mesh structure
 *
 */

structMesh* ParallelMeshCreate(){
    
    structMesh *mesh = (structMesh*)calloc(1, sizeof(structMesh));
    double xmin, xmax, ymin, ymax;
    int dim1, dim2;
    
    // start and end row number
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    int Rstart=0, Rend=0, Rnum;
    int p=0;
    int Rlocal=(np+nprocs-1)/nprocs;
    for (p=0; p<=procid; p++){
        Rend += Rlocal;
    }
    
    Rstart = Rend-Rlocal;
    if (Rend>np){
        Rend = np;
    }
    
    // Add two rows as ghost points.
    Rstart -= 1;
    Rend   += 1;
    Rnum = Rend-Rstart;
    
//    printf("procid=%d, Rstart=%d, Rlocal=%d, Rend=%d, Rnum=%d\n",
//           procid, Rstart, Rlocal, Rend, Rnum);
    
    // initialize the mesh
    xmin = 0.0; xmax = 2.0;
    ymin = 0.0; ymax = 2.0;
    
    mesh->ndim1 = Rnum; // No. of points along y coordinate
    mesh->ndim2 = np; // No. of points along x coordinate
    mesh->dx = (xmax - xmin)/(double)(np - 1.0);
    mesh->dy = (ymax - ymin)/(double)(np - 1.0);
    
    mesh->x = Matrix_create(mesh->ndim1, mesh->ndim2);
    mesh->y = Matrix_create(mesh->ndim1, mesh->ndim2);
    
    // initialize the condition
    for (dim1=0; dim1<Rnum; dim1++) {
        for (dim2=0; dim2<mesh->ndim2; dim2++) {
            int id = Rstart + dim1;
            mesh->x[dim1][dim2] = dim2*mesh->dx;
            mesh->y[dim1][dim2] = id*mesh->dy;
        }
    }
    
//    printf("procid=%d, Rstart=%d, Rend=%d, Rnum=%d, ymin=%f, ymax=%f\n",
//           procid, Rstart, Rend, Rnum, mesh->y[0][0], mesh->y[Rnum-1][0]);
    return mesh;
}

/**
 * @brief
 * Deallocate Memory for Mesh strucutre.
 */
void ParallelMeshFree(structMesh *mesh){
    Matrix_free(mesh->x);
    Matrix_free(mesh->y);
    return;
}

/**
 * @brief
 * Allocate Memory for concentration variable.
 */
double** ConcentrationCreate(structMesh *mesh, physics *phys){
    
    double **c = Matrix_create(mesh->ndim1, mesh->ndim2);
    double x0  = phys->x0;
    double y0  = phys->y0;
    double **x = mesh->x;
    double **y = mesh->y;
    double Dx  = phys->Dx;
    double Dy  = phys->Dy;
    
    int dim1, dim2;
    // initial condiation
    for (dim1=0; dim1< mesh->ndim1; dim1++) {
        for (dim2=0; dim2< mesh->ndim2; dim2++)
            c[dim1][dim2] = exp(-pow(x[dim1][dim2] - x0, 2.0)/Dx
                                - pow(y[dim1][dim2] - y0, 2.0)/Dy);
    }
    return c;
}

void ConcentrationFree(double **c){
    Matrix_free(c);
    return;
}

/**
 * @brief
 * Allocate Memory and Initizlize for Physics Structure.
 */
physics* PhysicsCreate(){
    
    physics *phys = (physics *)calloc(1, sizeof(physics));
    // initialize physics constant
    phys->u  = 1.0;
    phys->v  = 1.0;
    phys->Dx = 0.01;
    phys->Dy = 0.01;
    phys->x0 = 0.5;
    phys->y0 = 0.5;
    
    return phys;
}

/**
 * @brief
 * Deallocate Memory for physics structure.
 */
void PhysicsFree(physics *phys){
    
    free(phys);
    return;
}

/**
 * @brief
 * Exact Solution of Convection2d Problem.
 */
void ExactSol(structMesh *mesh, physics *phys, double time, double **c_exa){
    int dim1, dim2;
    double cx, cy;
    double **x = mesh->x;
    double **y = mesh->y;
    double u   = phys->u;
    double v   = phys->v;
    double x0  = phys->x0;
    double y0  = phys->y0;
    
    for (dim1=0; dim1<mesh->ndim1; dim1++) {
        for (dim2=0; dim2<mesh->ndim2; dim2++) {
            cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phys->Dx*(4.0*time+1) ) );
            cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phys->Dy*(4.0*time+1) ) );
            c_exa[dim1][dim2] = cx*cy/(4.0*time+1.0);
        }
    }
    return;
}

/**
 * @brief
 * Calculate the Norm Error of Numerical Solution at Spicific Time.
 */
void NormErr(structMesh *mesh, physics *phys,
             double time, double **c,
             double *L2, double *Linf){
    
    double **temp = Matrix_create(mesh->ndim1, mesh->ndim2);
    ExactSol(mesh, phys, time, temp);
    
    double err    = 0;
    double maxerr = 0.0;
    int    dim1, dim2;
    for (dim1=0; dim1<mesh->ndim1; dim1++) {
        for (dim2=0; dim2<mesh->ndim2; dim2++) {
            
            double dc = c[dim1][dim2]-temp[dim1][dim2];
            maxerr = maxf(fabs(dc), maxerr);
            err += dc*dc;
        }
    }
    *L2   = sqrt(err/mesh->ndim1/mesh->ndim2);
    *Linf = maxerr;
    
    Matrix_free(temp);
    return;
}
