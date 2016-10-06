//
//  ConvectionRHSUpwind2d.c
//  Project1
//
//  Created by li12242 on 16/10/5.
//  Copyright (c) 2016年 li12242. All rights reserved.
//

#include "ConvectionParallel2d.h"

/**
 * @brief
 * Calculate the Right Hand Side with Central FD Scheme.
 *
 * @param[in] mesh
 * @param[in] phys
 * @param[in] C
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * rhs | double ** |
 *
 */
void ConvectionRHS2d_central(structMesh *mesh, physics *phys,
                             double** C, double **rhs){
    
    int dim2, dim1;
    
    double u   = phys->u;
    double v   = phys->v;
    
    // rhs = -u * dc/dx
    for ( dim1=1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2 = 1; dim2<mesh->ndim2-1; dim2++) {
            double dxl = mesh->x[dim1  ][dim2  ] - mesh->x[dim1  ][dim2-1];
            double dxr = mesh->x[dim1  ][dim2+1] - mesh->x[dim1  ][dim2  ];
            
//            double c = C[dim1  ][dim2  ];
            double l = C[dim1  ][dim2-1];
            double r = C[dim1  ][dim2+1];
            // rhs = u*[(sign(u)+1)/2*(c[i][j] - c[i][j-1])/dx + (1-sign(u))/2*(c[i][j+1] - c[i][j])/dx ]
            
            rhs[dim1][dim2] = -u*(r-l)/(dxl+dxr);
        }
    }
    // rhs = -v * dc/dy + rhs
    for ( dim1=1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2 = 1; dim2<mesh->ndim2-1; dim2++) {
            double dyb = mesh->y[dim1  ][dim2  ] - mesh->y[dim1-1][dim2  ];
            double dyt = mesh->y[dim1+1][dim2  ] - mesh->y[dim1  ][dim2  ];
            
//            double c = C[dim1  ][dim2  ];
            double b = C[dim1-1][dim2  ];
            double t = C[dim1+1][dim2  ];
            
            
            rhs[dim1][dim2] = -v*(t-b)/(dyb+dyt);
        }
    }
    // rhs = Dx * d2c/dx2 + rhs
    for (dim1 = 1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2 = 1; dim2<mesh->ndim2-1; dim2++) {
            double dxl = mesh->x[dim1  ][dim2  ] - mesh->x[dim1  ][dim2-1];
            double dxr = mesh->x[dim1  ][dim2+1] - mesh->x[dim1  ][dim2  ];
            
            double c = C[dim1  ][dim2  ];
            double l = C[dim1  ][dim2-1];
            double r = C[dim1  ][dim2+1];
            
            rhs[dim1][dim2] += phys->Dx*( l+r-2*c )/(dxl*dxr);
        }
    }
    // rhs = Dy * d2c/dy2 + rhs
    for (dim1 = 1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2 = 1; dim2<mesh->ndim2-1; dim2++) {
            double dyb = mesh->y[dim1  ][dim2  ] - mesh->y[dim1-1][dim2  ];
            double dyt = mesh->y[dim1+1][dim2  ] - mesh->y[dim1  ][dim2  ];
            
            double c = C[dim1  ][dim2  ];
            double b = C[dim1-1][dim2  ];
            double t = C[dim1+1][dim2  ];
            
            rhs[dim1][dim2] += phys->Dx*( t+b-2*c )/(dyb*dyt);
        }
    }
    return;
}

/**
 * @brief
 * Calculate the Right Hand Side with Upwind FD Scheme.
 *
 * @param[in] mesh
 * @param[in] phys
 * @param[in] C
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * rhs | double ** |
 *
 */
void ConvectionRHS2d_upwind(structMesh *mesh, physics *phys,
                            double** C, double **rhs){
    
    int dim2, dim1;
    
    double u   = phys->u;
    double v   = phys->v;
    
    // rhs = -u * dc/dx
    for ( dim1=1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2 = 1; dim2<mesh->ndim2-1; dim2++) {
            double dxl = mesh->x[dim1  ][dim2  ] - mesh->x[dim1  ][dim2-1];
            double dxr = mesh->x[dim1  ][dim2+1] - mesh->x[dim1  ][dim2  ];
            
            double c = C[dim1  ][dim2  ];
            double l = C[dim1  ][dim2-1];
            double r = C[dim1  ][dim2+1];
            // rhs = u*[(sign(u)+1)/2*(c[i][j] - c[i][j-1])/dx + (1-sign(u))/2*(c[i][j+1] - c[i][j])/dx ]
            if (u>0) {
                rhs[dim1][dim2] = -u*(c-l)/dxl;
            }else{
                rhs[dim1][dim2] = -u*(r-c)/dxr;
            }
        }
    }
    // rhs = -v * dc/dy + rhs
    for ( dim1=1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2 = 1; dim2<mesh->ndim2-1; dim2++) {
            double dyb = mesh->y[dim1  ][dim2  ] - mesh->y[dim1-1][dim2  ];
            double dyt = mesh->y[dim1+1][dim2  ] - mesh->y[dim1  ][dim2  ];
            
            double c = C[dim1  ][dim2  ];
            double b = C[dim1-1][dim2  ];
            double t = C[dim1+1][dim2  ];
            
            if (v>0) {
                rhs[dim1][dim2] = -v*(c-b)/dyb;
            }else{
                rhs[dim1][dim2] = -v*(t-c)/dyt;
            }
        }
    }
    // rhs = Dx * d2c/dx2 + rhs
    for (dim1 = 1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2 = 1; dim2<mesh->ndim2-1; dim2++) {
            double dxl = mesh->x[dim1  ][dim2  ] - mesh->x[dim1  ][dim2-1];
            double dxr = mesh->x[dim1  ][dim2+1] - mesh->x[dim1  ][dim2  ];
            
            double c = C[dim1  ][dim2  ];
            double l = C[dim1  ][dim2-1];
            double r = C[dim1  ][dim2+1];
            
            rhs[dim1][dim2] += phys->Dx*( l+r-2*c )/(dxl*dxr);
        }
    }
    // rhs = Dy * d2c/dy2 + rhs
    for (dim1 = 1; dim1<mesh->ndim1-1; dim1++) {
        for (dim2 = 1; dim2<mesh->ndim2-1; dim2++) {
            double dyb = mesh->y[dim1  ][dim2  ] - mesh->y[dim1-1][dim2  ];
            double dyt = mesh->y[dim1+1][dim2  ] - mesh->y[dim1  ][dim2  ];
            
            double c = C[dim1  ][dim2  ];
            double b = C[dim1-1][dim2  ];
            double t = C[dim1+1][dim2  ];
            
            rhs[dim1][dim2] += phys->Dx*( t+b-2*c )/(dyb*dyt);
        }
    }
    return;
}