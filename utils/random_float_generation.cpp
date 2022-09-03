
#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"

 
// 产生[0,1]的浮点数，可以使用 rand()/double(RAND_MAX);
void random_double_generation(double* random, int num, double min, double max)
{

    srand((unsigned)time(NULL));
       
    for (int i = 0; i < num; i++)
		random[i] = (rand()/double(RAND_MAX) * (max - min) + min);
    
}


/*

int main (int argc, char const* argv[])
{
    double* M = new double[5];
    
	random_double_generation(M, 5, 1, 10);
	
	for(int j=0;j<5;j++)
		cout << M[j] << endl;
	
	return 0;
}

*/
