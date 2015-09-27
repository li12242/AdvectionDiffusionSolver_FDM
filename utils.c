//
// Created by li12242 on 15-9-14.
//
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

double *BuildVector(int Nrows){

    double *A = (double*) calloc(Nrows, sizeof(double));
    
    return A;
}

double *DestroyVector(double *v){
    free(v);
    return NULL;
}

int *BuildIntVector(int Nrows){

    int *A = (int*) calloc(Nrows, sizeof(int));

    return A;
}

int *DestroyIntVector(int *v){
    free(v);
    return NULL;
}


/* row major storage for a 2D matrix array */
int **BuildIntMatrix(int Nrows, int Ncols){
    int n;
    int **A = (int**) calloc(Nrows, sizeof(int*));

    A[0] = (int*) calloc(Nrows*Ncols, sizeof(int));

    for(n=1;n<Nrows;++n){
        A[n] = A[n-1]+ Ncols;
    }
    return A;
}

int **DestroyIntMatrix(int **A){
    free(A[0]);
    free(A);

    return NULL;
}


double **BuildMatrix(int Nrows, int Ncols){
    int n;
    double **A = (double**) calloc(Nrows, sizeof(double*));

    A[0] = (double*) calloc(Nrows*Ncols, sizeof(double));

    for(n=1;n<Nrows;++n){
        A[n] = A[n-1]+ Ncols;
    }

    return A;
}

double **DestroyMatrix(double **A){
    free(A[0]);
    free(A);

    return NULL;
}


void PrintMatrix(char *message, double **A, int Nrows, int Ncols){
    int n,m;

    printf("%s\n", message);
    for(n=0;n<Nrows;++n){
        for(m=0;m<Ncols;++m){
            printf(" %g ", A[n][m]);
        }
        printf(" \n");
    }
}


void SaveMatrix(char *filename, double **A, int Nrows, int Ncols){
    int n,m;

    FILE *fp = fopen(filename, "w");

    for(n=0;n<Nrows;++n){
        for(m=0;m<Ncols;++m){
            fprintf(fp, " %g ", A[n][m]);
        }
        fprintf(fp, " \n");
    }

    fclose(fp);
}


