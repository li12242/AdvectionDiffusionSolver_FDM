#include "ConvectionGPU2d.h"


//private func
void PrintTotalError(structMesh *mesh, physics *phys, float **c, float time);

int main(int argc, char **argv){
    
    physics    *phys = PhysicsCreate();
    structMesh *mesh = MeshCreate();
    float      **c   = ConcentrationCreate(mesh, phys);
    float      finalTime = 0.5f;
    
    PrintTotalError(mesh, phys, c, 0);
    ConvectionGPUSolve2d(mesh, phys, c, finalTime);
    // PrintTotalError(mesh, phys, c, 0);
    PrintTotalError(mesh, phys, c, finalTime);

    char filename[NAMELEN];
    snprintf(filename, NAMELEN, "result.txt");
//    printf("procid=%d, filename=%s\n", procid, filename);
    SaveMatrix(filename, c, mesh->ndim1, mesh->ndim2);
    snprintf(filename, NAMELEN, "y.txt");
    SaveMatrix(filename, mesh->y, mesh->ndim1, mesh->ndim2);
    snprintf(filename, NAMELEN, "x.txt");
    SaveMatrix(filename, mesh->x, mesh->ndim1, mesh->ndim2);
    
    free(phys);
    Matrix_free(c);
    MeshFree(mesh);
    return 0;
}

/* Print the Norm Error on screen */
void PrintTotalError(structMesh *mesh, physics *phys,
                     float **c, float time){
    float L2, Linf;
    
    NormErr(mesh, phys, time, c, &L2, &Linf);
    //printf("procid=%d, Local Norm Error, L2=%f, Linf=%f\n",procid,L2,Linf);
    printf("Total Norm Error, L2=%f, Linf=%f\n", L2,Linf);
    return;
}

/**
 * @brief
 * Allocate Memory for Mesh strucutre.
 *
 * @details
 * For multiply dimensional matrix, each dimensions represent 
 * the spacial coordinate x, y.
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh  | structMesh* | pointer to the mesh structure
 *
 */

structMesh* MeshCreate(){
    
    structMesh *mesh = (structMesh*)calloc(1, sizeof(structMesh));
    float xmin, xmax, ymin, ymax;
    int dim1, dim2;
    
    // initialize the mesh
    xmin = 0.0f; xmax = 2.0f;
    ymin = 0.0f; ymax = 2.0f;
    
    mesh->ndim1 = np; // No. of points along x coordinate
    mesh->ndim2 = np; // No. of points along y coordinate
    mesh->dx = (xmax - xmin)/(float)(np - 1.0f);
    mesh->dy = (ymax - ymin)/(float)(np - 1.0f);
    
    mesh->x = Matrix_create(mesh->ndim1, mesh->ndim2);
    mesh->y = Matrix_create(mesh->ndim1, mesh->ndim2);
    
    // initialize the condition
    for (dim1=0; dim1<mesh->ndim1; dim1++) {
        for (dim2=0; dim2<mesh->ndim2; dim2++) {
            mesh->x[dim1][dim2] = dim1*mesh->dx;
            mesh->y[dim1][dim2] = dim2*mesh->dy;
        }
    }
    return mesh;
}

/**
 * @brief
 * Deallocate Memory for Mesh strucutre.
 */
void MeshFree(structMesh *mesh){
    Matrix_free(mesh->x);
    Matrix_free(mesh->y);
    return;
}

/**
 * @brief
 * Allocate Memory for concentration variable.
 */
float** ConcentrationCreate(structMesh *mesh, physics *phys){
    
    float **c = Matrix_create(mesh->ndim1, mesh->ndim2);
    float x0  = phys->x0;
    float y0  = phys->y0;
    float **x = mesh->x;
    float **y = mesh->y;
    float Dx  = phys->Dx;
    float Dy  = phys->Dy;
    
    int dim1, dim2;
    // initial condiation
    for (dim1=0; dim1< mesh->ndim1; dim1++) {
        for (dim2=0; dim2< mesh->ndim2; dim2++)
            c[dim1][dim2] = exp(-pow(x[dim1][dim2] - x0, 2.0f)/Dx
                                - pow(y[dim1][dim2] - y0, 2.0f)/Dy);
    }
    return c;
}

void ConcentrationFree(float **c){
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
    phys->u  = 1.0f;
    phys->v  = 1.0f;
    phys->Dx = 0.01f;
    phys->Dy = 0.01f;
    phys->x0 = 0.5f;
    phys->y0 = 0.5f;

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
void ExactSol(structMesh *mesh, physics *phys, float time, float **c_exa){
    int dim1, dim2;
    float cx, cy;
    float **x = mesh->x;
    float **y = mesh->y;
    float u   = phys->u;
    float v   = phys->v;
    float x0  = phys->x0;
    float y0  = phys->y0;
    
    for (dim1=0; dim1<mesh->ndim1; dim1++) {
        for (dim2=0; dim2<mesh->ndim2; dim2++) {
            cx = exp( -pow( x[dim1][dim2]-x0-u*time, 2.0f )/( phys->Dx*(4.0f*time+1) ) );
            cy = exp( -pow( y[dim1][dim2]-y0-v*time, 2.0f )/( phys->Dy*(4.0f*time+1) ) );
            c_exa[dim1][dim2] = cx*cy/(4.0f*time+1.0);
        }
    }
    return;
}

/**
 * @brief
 * Calculate the Norm Error of Numerical Solution at Spicific Time.
 */
void NormErr(structMesh *mesh, physics *phys,
             float time, float **c,
             float *L2, float *Linf){
    
    float **temp = Matrix_create(mesh->ndim1, mesh->ndim2);
    ExactSol(mesh, phys, time, temp);

    float err    = 0.0f;
    float maxerr = 0.0f;
    int   dim1,dim2;
    
    for (dim1=0; dim1<mesh->ndim1; dim1++) {
        for (dim2=0; dim2<mesh->ndim2; dim2++) {
            
            float dc = c[dim1][dim2]-temp[dim1][dim2];
            maxerr = maxf(fabs(dc), maxerr);
            err += dc*dc;
        }
    }
    *L2   = sqrt(err/mesh->ndim1/mesh->ndim2);
    *Linf = maxerr;
    
    Matrix_free(temp);
    return;
}
