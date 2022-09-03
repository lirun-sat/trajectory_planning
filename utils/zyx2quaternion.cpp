#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"


void zyx2quaternion(double psi, double theta, double phi, double* quaternion)
{
    
	/*
	 欧拉角转四元素，其中 欧拉角 用 动系相继旋转 Z-Y-X法表示
     psi: 表示 Z 轴转动量
     theta: 表示 Y 轴转动量
     phi: 表示 X 轴转动量
     return: 四元素，q0标量，q1,2,3向量
	*/
    
	quaternion[0] = cos(phi/2) * cos(theta/2) * cos(psi/2) + sin(phi/2) * sin(theta/2) * sin(psi/2);
	
	quaternion[1] = sin(phi/2) * cos(theta/2) * cos(psi/2) - cos(phi/2) * sin(theta/2) * sin(psi/2);
	
	quaternion[2] = cos(phi/2) * sin(theta/2) * cos(psi/2) + sin(phi/2) * cos(theta/2) * sin(psi/2);
	
	quaternion[3] = cos(phi/2) * cos(theta/2) * sin(psi/2) - sin(phi/2) * sin(theta/2) * cos(psi/2);
    
    
}


/*
int main()
{
	double z = 45 * 3.14 / 180;
	double y = -50 * 3.14 / 180;
	double x = -30 * 3.14 / 180;
	
    double* M = new double[4];
    zyx2quaternion(z, y, x, M);
    
	for(int i = 0; i < (4); i++)
		cout << M[i] << endl;
    
    delete[] M;

    return 0;
}
*/













