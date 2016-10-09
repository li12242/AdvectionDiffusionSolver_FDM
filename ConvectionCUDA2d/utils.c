//
// Created by li12242 on 15-9-14.
//
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

float *Vector_create(int Nrows){
    float *A = (float*) calloc(Nrows, sizeof(float));
    return A;
}

float *Vector_free(float *v){
    free(v);
    return NULL;
}

int *IntVector_create(int Nrows){
    int *A = (int*) calloc(Nrows, sizeof(int));
    return A;
}

int *IntVector_free(int *v){
    free(v);
    return NULL;
}


/* row major storage for a 2D matrix array */
int **IntMatrix_create(int Nrows, int Ncols){
    int n;
    int **A = (int**) calloc(Nrows, sizeof(int*));

    A[0] = (int*) calloc(Nrows*Ncols, sizeof(int));

    for(n=1;n<Nrows;++n){
        A[n] = A[n-1]+ Ncols;
    }
    return A;
}

int **IntMatrix_free(int **A){
    free(A[0]);
    free(A);

    return NULL;
}


float **Matrix_create(int Nrows, int Ncols){
    int n;
    float **A = (float**) calloc(Nrows, sizeof(float*));

    A[0] = (float*) calloc(Nrows*Ncols, sizeof(float));

    for(n=1;n<Nrows;++n){
        A[n] = A[n-1]+ Ncols;
    }

    return A;
}

float **Matrix_free(float **A){
    free(A[0]);
    free(A);

    return NULL;
}


void PrintMatrix(char *message, float **A, int Nrows, int Ncols){
    int n,m;

    printf("%s\n", message);
    for(n=0;n<Nrows;++n){
        for(m=0;m<Ncols;++m){
            printf(" %g ", A[n][m]);
        }
        printf(" \n");
    }
}


void SaveMatrix(char *filename, float **A, int Nrows, int Ncols){
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


