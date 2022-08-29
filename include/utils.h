#ifndef CALC_UTILS
#define CALC_UTILS


#include<stdio.h>
#include<cmath>
#include<iostream>
#include<ctime>
#include<cstring>
#include<cstdlib>

using namespace std;



void random_double_generation(double* random, int num, double min, double max);

void MatrixMulti(double* A, double* B, double* AB);

void MatrixMulti_(int row_1, int col_1, int col_2, double* A, double* B, double* AB);

void zyx2quaternion(double psi, double theta, double phi, double* quaternion);

void quaternion2dc(double q_0, double q_1, double q_2, double q_3, double* dc);

void rpy2dc(double roll, double pitch, double yaw, double* direction_cosine);

void dc2euler(double* R, double* euler);

void euler_ZYX2dc(double alpha, double beta, double gamma, double* A);

void quaternion2euler(double w, double x, double y, double z, double* euler);

void MatrixMultiTranspose(double* A, double* B, double* AB);

void MatrixTranspose_(int row, int col, double* M, double* M_T);

void MatrixTranspose(double* M);

void MatrixAdd_(int row, int col, double* A, double* B, double* AB);

void MatrixSub_(int row, int col, double* A, double* B, double* AB);

void ScaleMatrix_(int row, int col, double scale, double* A, double* scale_A);

void SetZeroMatrix_(int row, int col, double* A);

void MatrixExtract_(int m, int n, int i_s, int i_f, int j_s, int j_f, double *A, double *ans);

void MatrixCopy_(int m, int n, double *a, double *u);

void MatrixDiagExpand(double* A, int row, int col, double* A_expand);

void cross(double* V, double* M);

void cxyz(double theta, int x, int y, int z, double* direction_cosines);

void LUP_Descomposition(double* A, double* L, double* U, int* P, int n);

void LUP_Solve(double* L, double* U, int* P, double* b, int n, double* inv_A_each);

int getNext(int i, int m, int n);

int getPre(int i, int m, int n);

void movedata(double *mtx, int i, int m, int n);

void transpose(double *mtx, int m, int n);

void LUP_solve_inverse(double* A, int n, double* inv_A);

double calc_determinantOfMatrix(double* a, int num, int row, int col, int n);

void calc_Euler_dev2angvl(double* Euler_angle, double* Euler_dev, double* angvl);


#endif 
