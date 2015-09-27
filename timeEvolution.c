//
// Created by li12242 on 15-9-15.
//

#include "timeEvolution.h"

double calTimeStep(structMesh *mesh, physics *phy, double);
int calRHS_upwind(structMesh *, physics *, double**, double **);
int updateBC(structMesh *mesh, physics *phy, double time, double ** concentration);

/** \fn int timeEvolution(structMesh *mesh, physics *phy, double*** concentration)
 \brief use Euler forward time evolution
 \param   mesh           mesh structure
 \param   phy               physical constant
 \param   concentration    variable
 \warning
 */

int timeEvolution(structMesh *mesh, physics *phy, double** concentration){
    double time, CFL, deltaT, finalTime;
    double **rhs = BuildMatrix(mesh->ny, mesh->nx);
    const double *temp = rhs[0];
    
    finalTime = 0.5;
    time = 0.0;
    CFL = 0.1;
    
    temp = rhs[0];
    // cal delta T
    deltaT = calTimeStep(mesh, phy, CFL);
    
    while (time < finalTime) {
        
        // increase time
        if (time+deltaT > finalTime) {
            time = finalTime;
        }else{
            time += deltaT;
        }
        
        updateBC(mesh, phy, time, concentration);
        calRHS_upwind(mesh, phy, concentration, rhs);
        
        // concentration = deltaT * rhs + concentration
        cblas_daxpy( mesh->nx*mesh->ny, deltaT, *rhs, 1, *concentration, 1);
        
    }
    
    DestroyMatrix(rhs);
    return 0;
}

/** \fn  double calTimeStep(structMesh *mesh, physics *phy, double CFL)
 \brief  calculate mininum delta time
 \param  mesh
 \param  phy
 \param  CFL
 \warning
 */

double calTimeStep(structMesh *mesh, physics *phy, double CFL){
    double deltaT;
    deltaT = mesh->dx/phy->u;
    deltaT = minf(deltaT, mesh->dy/phy->v);
    deltaT = minf(deltaT, pow(mesh->dx, 2)/phy->Dx);
    deltaT = minf(deltaT, pow(mesh->dy, 2)/phy->Dy);
    return deltaT*CFL;
}

/** \fn  int calRHS(structMesh *mesh, physics *phy, double ***rhs)
 \brief  calcation Right Hands Side
 \param  mesh
 \param  phy
 \param  rhs
 \warning
 */

int calRHS_upwind(structMesh *mesh, physics *phy, double** concentration, double **rhs){
    
    int icol, irow;
    int sign;
    double **temp = BuildMatrix(mesh->ny, mesh->nx); /// \var temp description
    // rhs = -u * dc/dx
    sign = signf(phy->u);
    for ( irow=1; irow<mesh->ny-1; irow++) {
        for (icol = 1; icol<mesh->nx-1; icol++) {
            // rhs = u*[(sign(u)+1)/2*(c[i][j] - c[i][j-1])/dx + (1-sign(u))/2*(c[i][j+1] - c[i][j])/dx ]
            rhs[irow][icol] = -phy->u*( (sign/2+0.5)*(concentration[irow][icol] - concentration[irow][icol-1])/mesh->dx
                                       + (0.5-sign/2)*(concentration[irow][icol+1] - concentration[irow][icol])/mesh->dx );
        }
    }
    // rhs = -v * dc/dy + rhs
    sign = signf(phy->v);
    for ( irow=1; irow<mesh->ny-1; irow++) {
        for (icol = 1; icol<mesh->nx-1; icol++) {
            rhs[irow][icol] = rhs[irow][icol] -phy->v*( (sign/2+0.5)*(concentration[irow][icol] - concentration[irow-1][icol])/mesh->dy
                                          + (0.5-sign/2)*(concentration[irow+1][icol] - concentration[irow][icol])/mesh->dy );
        }
    }
    // rhs = Dx * d2c/dx2 + rhs
    for (irow = 1; irow<mesh->ny-1; irow++) {
        for (icol = 1; icol<mesh->nx-1; icol++) {
            rhs[irow][icol] = rhs[irow][icol] +
            phy->Dx*( concentration[irow][icol-1] + concentration[irow][icol+1] - 2*concentration[irow][icol] )/(mesh->dx*mesh->dx);
        }
    }
    // rhs = Dy * d2c/dy2 + rhs
    for (irow = 1; irow<mesh->ny-1; irow++) {
        for (icol = 1; icol<mesh->nx-1; icol++) {
            rhs[irow][icol] = rhs[irow][icol] +
            phy->Dx*( concentration[irow+1][icol] + concentration[irow-1][icol] - 2*concentration[irow][icol] )/(mesh->dy*mesh->dy);
        }
    }
    DestroyMatrix(temp);
    return 0;
}

/** \fn  int updateBC(structMesh *mesh, physics *phy, double ** concentration)
 \brief  update 4 boundaries
 \param  mesh mesh structure
 \param  phy  physical constant
 \param  concentration  variable
 \warning
 */
int updateBC(structMesh *mesh, physics *phy, double time, double ** concentration){
    int irow, icol;
    
    for (irow = 0; irow<mesh->ny; irow++) {
        for (icol = 0; icol<mesh->nx; icol+=(mesh->nx-1) ) {
            concentration[irow][icol] = exp(-pow((mesh->x[icol]-0.5-phy->u*time), 2.0)/(phy->Dx*(4*time+1)))
            *exp( -pow(mesh->y[irow] - 0.5 -phy->v*time, 2.0)/(phy->Dy*(4*time+1)) )/(4.0*time+1.0);
        }
    }
    
    for (irow=0; irow<mesh->ny; irow+=(mesh->ny-1)) {
        for (icol = 0; icol<mesh->nx; icol++) {
            concentration[irow][icol] = exp(-pow((mesh->x[icol]-0.5-phy->u*time), 2.0)/(phy->Dx*(4*time+1)))
            *exp( -pow(mesh->y[irow] - 0.5 -phy->v*time, 2.0)/(phy->Dy*(4*time+1)) )/(4.0*time+1.0);
        }
    }
    return 0;
}



