#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void euler_ZYX2dc(double alpha, double beta, double gamma, double* A)
{
	// A 是一个 一维数组（3 * 3），用这个一维数组来表示二维数组，元素用 i * N + j 来查找
    // Z-Y-X 欧拉角[alpha, beta, gamma]定义为：frame_j分别绕其 Z 轴、旋转后的 Y 轴、再旋转后的 X 轴旋转alpha, beta, gamma角后, 与frame_i重合
	A[0 * 3 + 0] = cos(alpha) * cos(beta);
	A[0 * 3 + 1] = cos(alpha) * sin(beta) * sin(gamma) - sin(alpha) * cos(gamma);
    A[0 * 3 + 2] = cos(alpha) * sin(beta) * cos(gamma) + sin(alpha) * sin(gamma);

    A[1 * 3 + 0] = sin(alpha) * cos(beta);
    A[1 * 3 + 1] = sin(alpha) * sin(beta) * sin(gamma) + cos(alpha) * cos(gamma);
    A[1 * 3 + 2] = sin(alpha) * sin(beta) * cos(gamma) - cos(alpha) * sin(gamma);

    A[2 * 3 + 0] = -sin(beta);
    A[2 * 3 + 1] = cos(beta) * sin(gamma);
    A[2 * 3 + 2] = cos(beta) * cos(gamma);
	
}


/*
int main()
{
    double alpha;
    double beta;
	double gamma;
	
	alpha = -150 * PI / 180;
	beta = 150 * PI / 180;
	gamma = -150 * PI / 180;

    double* M = new double[9];
	
    euler_ZYX2dc(alpha, beta, gamma, M);

    cout << M[0] << " " << M[1] << " " << M[2] << endl;
    cout << M[3] << " " << M[4] << " " << M[5] << endl;
	cout << M[6] << " " << M[7] << " " << M[8] << endl;
    delete[] M;

    return 0;
}

*/













