#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "math.h"
#include "timeEvolution.h"

#define np 101 // number of points along single axis

int InitConst(structMesh *, physics *);
//double** InitVar(structMesh*, physics*, double **);
int InitVar(structMesh*, physics*, double ***);
int Finalization(structMesh *, double **);

int main(){

    structMesh *mesh = (structMesh *)calloc(1, sizeof(structMesh)); // allocate memory
    physics *phy = (physics *)calloc(1, sizeof(physics));
    double **concentration = NULL;
    
    InitConst(mesh, phy);
    
    concentration = BuildMatrix(np, np);
    
//    concentration = InitVar(mesh, phy, concentration);
    InitVar(mesh, phy, &concentration);
    
    PrintMatrix("con1", concentration, mesh->ny, mesh->nx);
    
    timeEvolution(mesh, phy, concentration);
    
    PrintMatrix("con2", concentration, mesh->ny, mesh->nx);

    Finalization(mesh, concentration);
    return 0;
}


int InitConst(structMesh *mesh, physics *phy){
    double xmin, xmax, ymin, ymax;
    int icol, irow;
    // initialize the mesh
    
    xmin = 0.0; xmax = 2.0;
    ymin = 0.0; ymax = 2.0;
    
    mesh->nx = np;
    mesh->ny = np;
    mesh->dx = (xmax - xmin)/(double)(mesh->nx - 1);
    mesh->dy = (ymax - ymin)/(double)(mesh->ny - 1);
    
    mesh->x = BuildVector(mesh->nx);
    mesh->y = BuildVector(mesh->ny);
    
    for (icol=0; icol< mesh->nx-1; icol++) {
        mesh->x[icol+1] = mesh->x[icol] + mesh->dx;
    }
    
    for (irow=0; irow< mesh->ny-1; irow++) {
        mesh->y[irow+1] = mesh->y[irow] + mesh->dy;
    }
    
    // initialize physics constant
    phy->u = 1.0;
    phy->v = 1.0;
    phy->Dx = 0.01;
    phy->Dy = 0.01;

    return 0;
}


//double** InitVar(structMesh *mesh, physics *phy, double ***concentration){
int InitVar(structMesh *mesh, physics *phy, double ***concentration){
    int irow, icol;

//  set the initial condiation
    
    
    for (irow=0; irow< mesh->ny; irow++) {
        for (icol=0; icol< mesh->nx; icol++) {
            (*concentration)[irow][icol] = exp(-pow(mesh->x[icol] - 0.5, 2.0)/phy->Dx
                                            - pow(mesh->y[irow] - 0.5, 2.0)/phy->Dy);
        }
    }

//    return concentration;
    return 0;
}


int Finalization(structMesh *mesh, double **concentration){
    DestroyVector(mesh->x);
    DestroyVector(mesh->y);
    DestroyMatrix(concentration);
    return 0;
}