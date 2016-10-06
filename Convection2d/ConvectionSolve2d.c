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
    
    time      = 0.0;
    CFL       = 0.1;
    // cal dt
    dt = getTimeInterval(mesh, phy, CFL);
    
    while (time < finalTime) {
        
        // increase time
        if (time+dt > finalTime) {
            dt   = finalTime - time;
            time = finalTime;
        }else{
            time += dt;
        }
        
        ConvectionBC2d(mesh, phy, time, c);
//        ConvectionRHS2d_central(mesh, phy, c, rhs);
        ConvectionRHS2d_upwind(mesh, phy, c, rhs);
        
        // c = deltaT * rhs + c
        int dim1, dim2;
        for (dim1=0; dim1<mesh->ndim1; dim1++) {
            for (dim2=0; dim2<mesh->ndim2; dim2++) {
                c[dim1][dim2] += dt*rhs[dim1][dim2];
            }
        }
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
    double dx = mesh->dx;
    double dy = mesh->dy;
    dt = dx/phy->u;
    dt = minf(dt, dy/phy->v);
    dt = minf(dt, dx*dx/phy->Dx);
    dt = minf(dt, dy*dy/phy->Dy);
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
        // left side
        dim2 = 0;
        cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phy->Dx*(4.0*time+1) ) );
        cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phy->Dy*(4.0*time+1) ) );
        c[dim1][dim2] = cx*cy/(4.0*time+1.0);
        // right side
        dim2 = mesh->ndim2-1;
        cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phy->Dx*(4.0*time+1) ) );
        cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phy->Dy*(4.0*time+1) ) );
        c[dim1][dim2] = cx*cy/(4.0*time+1.0);
    }
    
    for (dim2=0; dim2<mesh->ndim2; dim2++) {
        // bottom
        dim1 = 0;
        cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phy->Dx*(4.0*time+1) ) );
        cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phy->Dy*(4.0*time+1) ) );
        c[dim1][dim2] = cx*cy/(4.0*time+1.0);
        // top
        dim1 = mesh->ndim1-1;
        cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0 )/( phy->Dx*(4.0*time+1) ) );
        cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0 )/( phy->Dy*(4.0*time+1) ) );
        c[dim1][dim2] = cx*cy/(4.0*time+1.0);
    }
    return;
}





