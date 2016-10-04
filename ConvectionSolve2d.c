//
// Created by li12242 on 15-9-15.
//

#include "Convection2d.h"

/* private function */
double getTimeInterval(structMesh *mesh, physics *phys, double CFL);
void   ConvectionBC2d (structMesh *mesh, physics *phys, double time, double **c);

/**
 * @brief 
 * Use Euler Advection Scheme to Solve the PDE.
 */
int ConvectionSolve2d(structMesh *mesh, physics *phy, double** c, double finalTime){
    double time, CFL, dt;
    double **rhs = Matrix_create(mesh->ndim1, mesh->ndim2);
    const double *temp = rhs[0];
    
    finalTime = 0.5;
    time      = 0.0;
    CFL       = 0.1;
    
    temp = rhs[0];
    // cal dt
    dt = getTimeInterval(mesh, phy, CFL);
    
    while (time < finalTime) {
        
        // increase time
        if (time+dt > finalTime) {
            time = finalTime;
        }else{
            time += dt;
        }
        
        ConvectionBC2d(mesh, phy, time, c);
//        ConvectionRHS2d_central(mesh, phy, c, rhs);
        ConvectionRHS2d_upwind(mesh, phy, c, rhs);
        
        // c = deltaT * rhs + c
        cblas_daxpy( mesh->ndim2*mesh->ndim1, dt, *rhs, 1, *c, 1);
    }
    
    Matrix_free(rhs);
    return 0;
}

/**
 * @brief
 * Calculate the Time Interval.
 *
 */
double getTimeInterval(structMesh *mesh, physics *phy, double CFL){
    double dt;
    dt = mesh->dx/phy->u;
    dt = minf(dt, mesh->dy/phy->v);
    dt = minf(dt, pow(mesh->dx, 2)/phy->Dx);
    dt = minf(dt, pow(mesh->dy, 2)/phy->Dy);
    return dt*CFL;
}



/**
 * @brief
 * Update Boundary Values through Analytic Solutions.
 */
void ConvectionBC2d(structMesh *mesh, physics *phy,
                    double time, double **c){
    
    int dim1, dim2;
    double cx, cy;
    double **x = mesh->x;
    double **y = mesh->y;
    double u   = phy->u;
    double v   = phy->v;
    double x0  = phy->x0;
    double y0  = phy->y0;
    
    for (dim1 = 0; dim1<mesh->ndim1; dim1++) {
        dim2 = 0;
        cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phy->Dx*(4.0*time+1) ) );
        cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phy->Dy*(4.0*time+1) ) );
        c[dim1][dim2] = cx*cy/(4.0*time+1.0);
        
        dim2 = mesh->ndim2-1;
        cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phy->Dx*(4.0*time+1) ) );
        cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phy->Dy*(4.0*time+1) ) );
        c[dim1][dim2] = cx*cy/(4.0*time+1.0);
    }
    
    for (dim2=0; dim2<mesh->ndim2; dim2++) {
        dim1 = 0;
        cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phy->Dx*(4.0*time+1) ) );
        cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phy->Dy*(4.0*time+1) ) );
        c[dim1][dim2] = cx*cy/(4.0*time+1.0);
        
        dim1 = mesh->ndim1-1;
        cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phy->Dx*(4.0*time+1) ) );
        cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phy->Dy*(4.0*time+1) ) );
        c[dim1][dim2] = cx*cy/(4.0*time+1.0);
    }
    return;
}





