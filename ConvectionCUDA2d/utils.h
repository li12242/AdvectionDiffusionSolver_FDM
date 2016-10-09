//
// Created by li12242 on 15-9-14.
//

#ifndef PROGRAM_UTILS_H
#define PROGRAM_UTILS_H

#define signf(a)   ( (a>0)?1:-1 )
#define minf(a, b)  ( (a>b)?b:a )
#define maxf(a, b)  ( (a>b)?a:b )

float *Vector_create   (int Nrows);
float *Vector_free (float *v);
int    *IntVector_create(int Nrows);
int    *IntVector_free (int *v);
int    **IntMatrix_create  (int Nrows, int Ncols);
int    **IntMatrix_free(int **A);
float **Matrix_create     (int Nrows, int Ncols);
float **Matrix_free   (float **A);

void   PrintMatrix(char *message, float **A, int Nrows, int Ncols);
void   SaveMatrix (char *filename, float **A, int Nrows, int Ncols);

#endif //PROGRAM_UTILS_H