//
// Created by li12242 on 15-9-15.
//

#include "ConvectionParallel2d.h"

/* private function */
double getTimeInterval(structMesh *mesh, physics *phys, double CFL);
void   ConvectionBC2d (structMesh *mesh, physics *phys, double time, double **c);
void   TimeEvolution(structMesh *mesh, physics *phys,
                     double **c, double time, double dt, double **rhs);
void   ConvectionMPISendRecv2d(structMesh *mesh, double **c);
/**
 * @brief 
 * Use Euler Advection Scheme to Solve the PDE.
 */
int ConvectionSolve2d(structMesh *mesh, physics *phy, double** c, double finalTime){
    double time, CFL, dt;
    double **rhs = Matrix_create(mesh->ndim1, mesh->ndim2);
    
    time      = 0.0;
    CFL       = 0.1;
    
    // calculate time interval - dt
    dt = getTimeInterval(mesh, phy, CFL);
    
    while (time < finalTime) {
        
        // increase local time
        if (time+dt > finalTime) {
            dt   = finalTime - time;
            time = finalTime;
        }else{
            time += dt;
        }
        
        TimeEvolution(mesh, phy, c, time, dt, rhs);
    }
    
    Matrix_free(rhs);
    return 0;
}

/**
 * @brief 
 * Euler Time Advection Scheme.
 */
void TimeEvolution(structMesh *mesh, physics *phys,
                   double **c, double time, double dt, double **rhs){
    
    ConvectionBC2d(mesh, phys, time, c);
    ConvectionMPISendRecv2d(mesh, c);
    ConvectionRHS2d_upwind(mesh, phys, c, rhs);
    
    // c += dt * rhs
    int dim1, dim2;
    for (dim1=1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2=1; dim2<mesh->ndim2-1; dim2++) {
            c[dim1][dim2] += dt*rhs[dim1][dim2];
        }
    }
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
    return;
}

void ConvectionMPISendRecv2d(structMesh *mesh, double **c){
    double *CTopIn  = Vector_create(mesh->ndim2);
    double *CBotIn  = Vector_create(mesh->ndim2);
    double *CTopOut = Vector_create(mesh->ndim2);
    double *CBotOut = Vector_create(mesh->ndim2);
    
    int dim1, dim2;
    for (dim2=0; dim2<mesh->ndim2; dim2++) {
        dim1=1;
        CTopOut[dim2] = c[dim1][dim2];
        dim1=mesh->ndim1-2;
        CBotOut[dim2] = c[dim1][dim2];
    }
    
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    MPI_Request sendRequest[2], recvRequest[2];
    
    int topProcid=((procid-1)<0         )?0:(procid-1);
    int BotProcid=((procid+1)>(nprocs-1))?(nprocs-1):(procid+1);
    
//    printf("procid=%d, top procid=%d, bot procid=%d\n", procid, topProcid, BotProcid);
    
    // use destination process id as tags.
    MPI_Isend(CTopOut, mesh->ndim2, MPI_DOUBLE, topProcid, topProcid+665, MPI_COMM_WORLD, sendRequest);
    MPI_Isend(CBotOut, mesh->ndim2, MPI_DOUBLE, BotProcid, BotProcid, MPI_COMM_WORLD, sendRequest+1);
    
    // receive the message by own process id.
    int topTags=procid;
    int botTags=procid+665;
    if (procid==0) {
        topTags=procid+665; // for top process
    }else if(procid==nprocs-1){
        botTags=procid; // for bottom process
    }
    
    MPI_Irecv(CTopIn, mesh->ndim2, MPI_DOUBLE, topProcid, topTags, MPI_COMM_WORLD, recvRequest);
    MPI_Irecv(CBotIn, mesh->ndim2, MPI_DOUBLE, BotProcid, botTags, MPI_COMM_WORLD, recvRequest+1);
    
    MPI_Status instatus[2];
    MPI_Status outstatus[2];
    
    MPI_Waitall(2, recvRequest, instatus);
    MPI_Waitall(2, sendRequest, outstatus);
    
    for (dim2=0; dim2<mesh->ndim2; dim2++) {
        dim1=0;
        c[dim1][dim2] = CTopIn[dim2];
        dim1=mesh->ndim1-1;
        c[dim1][dim2] = CBotIn[dim2];
    }
    
    Vector_free(CTopIn);
    Vector_free(CTopOut);
    Vector_free(CBotIn);
    Vector_free(CBotOut);
    
    return;
}






