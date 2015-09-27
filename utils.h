//
// Created by li12242 on 15-9-14.
//

#ifndef PROGRAM_UTILS_H
#define PROGRAM_UTILS_H

double *BuildVector(int Nrows);

double *DestroyVector(double *v);

int *BuildIntVector(int Nrows);

int *DestroyIntVector(int *v);

int **BuildIntMatrix(int Nrows, int Ncols);

int **DestroyIntMatrix(int **A);

double **BuildMatrix(int Nrows, int Ncols);

double **DestroyMatrix(double **A);

void PrintMatrix(char *message, double **A, int Nrows, int Ncols);

void SaveMatrix(char *filename, double **A, int Nrows, int Ncols);

typedef struct{
    int nx;     // num of points along x axis
    int ny;     // num of points along y axis
    double dx;  // distance between x points
    double dy;  // distance between y points
    double *x;  // x coordinate
    double *y;  // y coordinate
} structMesh;

typedef struct{
    double u;
    double v;
    double Dx;
    double Dy;
} physics;

#define signf(a)  ( (a>0)?1:-1 )
#define minf(a, b)  ( (a>b)?b:a )

#endif //PROGRAM_UTILS_H