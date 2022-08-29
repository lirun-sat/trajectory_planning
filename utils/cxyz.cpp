#define _USE_MATH_DEFINES 
#include "utils.h"


void cxyz(double theta, int x, int y, int z, double* direction_cosines)
{
    if (x == 1)
	{
		direction_cosines[0] = 1; direction_cosines[1] = 0;          direction_cosines[2] = 0;
		
		direction_cosines[3] = 0; direction_cosines[4] = cos(theta); direction_cosines[5] = -sin(theta);
		
		direction_cosines[6] = 0; direction_cosines[7] = sin(theta); direction_cosines[8] = cos(theta);
		
	}
	else if (y == 1)
	{
		direction_cosines[0] = cos(theta);  direction_cosines[1] = 0; direction_cosines[2] = sin(theta);
		
		direction_cosines[3] = 0;           direction_cosines[4] = 1; direction_cosines[5] = 0;
		
		direction_cosines[6] = -sin(theta); direction_cosines[7] = 0; direction_cosines[8] = cos(theta);
		
	}
	else if (z == 1)
	{
		direction_cosines[0] = cos(theta); direction_cosines[1] = -sin(theta); direction_cosines[2] = 0;
		
		direction_cosines[3] = sin(theta); direction_cosines[4] = cos(theta);  direction_cosines[5] = 0;
		
		direction_cosines[6] = 0;          direction_cosines[7] = 0;           direction_cosines[8] = 1;
		
	}
    
}
/*
int main()
{
	double theta = 30 * 3.14 / 180;
	
	int x = 1;
	
    double* M = new double[3 * 3];
    cxyz(theta, x, 0, 0, M);
    
	for(int i = 0; i < (3*3); i++)
		cout << M[i] << endl;
    
    delete[] M;

    return 0;
}

*/






















