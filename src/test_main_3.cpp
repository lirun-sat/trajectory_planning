#include "WOAAlgorithm.h"
#include "PSOAlgorithm.h"

using namespace std;


// 采用 经典 函数 测试 woa_tlbo 算法 
//待优化的目标函数
double objFunction_PSO(PSO_Particle& particle);
double objFunction_WOA(WOA_Whale& whale);
void sort_max2min_main(double* , int , int* );
double avg(double* parameter, int n);
double stddev(double* parameter, int n);


int main(void)
{    
    double* minPos = new double[30];
    double* maxPos = new double[30];
    for(int i = 0; i < 30; i++)
    {
    	minPos[i] = -500;
    	maxPos[i] = 500;
    }
    int N = 30;
    int dimension = N;
    int pso_woa_count = 30;
	int Iter_Max = 800;
	int num_calc = 10;
	double inertGuideCoe = 0.7;
	double globalGuideCoe = 1.4961;
	double localGuideCoe = 1.4961;
	double disturbanceRate = 0;  // 0.2;
	double disturbanceVelocityCoe = 0;  // 0.05;
    
    double* tempt = new double[num_calc];  
	double* pso_woa_fitness = new double[pso_woa_count];
	int* pso_woa_index = new int[pso_woa_count];
	double** pso_woa_position = new double*[pso_woa_count];
	for( int i = 0; i < pso_woa_count; i++)
	{
		pso_woa_position[i] = new double[dimension];
	}
    
    double* maxSpeed = new double[N];
	double* minSpeed = new double[N];
    maxSpeed[0] = 1.5;
	maxSpeed[1] = 0.628;
	maxSpeed[2] = 1.5;
	maxSpeed[3] = 1.5;
	maxSpeed[4] = 1;
	maxSpeed[5] = 0.8;
	maxSpeed[6] = 1.5;

	minSpeed[0] = -1.5;
	minSpeed[1] = -0.628;
	minSpeed[2] = -1.5;
	minSpeed[3] = -1.5;
	minSpeed[4] = -1.5;
	minSpeed[5] = -0.8;
	minSpeed[6] = -1.5;
	

    //构建算法
    WOA_Algorithm woa(objFunction_WOA, minPos, maxPos, dimension, pso_woa_count/2);
	
	PSO_Algorithm pso(objFunction_PSO, minPos, maxPos, dimension, pso_woa_count / 2, maxSpeed, minSpeed, inertGuideCoe, globalGuideCoe, localGuideCoe);
	
	// *********************************************************************************************************
	// cout << "构建成功" << endl;
	
	for(int kk = 0; kk < num_calc; kk++)
	{
		woa.randomlyInitial();

		pso.randomlyInitial();
		
	    for(int Iter = 0; Iter < Iter_Max; Iter++)
	    {
	    	for(int i = 0; i < pso_woa_count; i++)
	    	{
	    		if(i < pso_woa_count/2)
	    		{
	    			pso_woa_fitness[i] = woa._whaleSet[i]._fitness;
    
	    			for(int j = 0; j < dimension; j++)
	    			{
	    				pso_woa_position[i][j] = woa._whaleSet[i]._position[j];
	    			}
 			
	    		}
	    		
	    		else
	    		{
	    			pso_woa_fitness[i] = pso._particleSet[i - pso_woa_count/2]._fitness;
	    		
	    			for(int j = 0; j < dimension; j++)
	    			{
	    				pso_woa_position[i][j] = pso._particleSet[i - pso_woa_count/2]._position[j];
	    			}

	    		}
                pso_woa_index[i] = i;
	    	 
	    	}            	
	    	    
	    	sort_max2min_main(pso_woa_fitness, pso_woa_count, pso_woa_index);
	    
	    	for(int i = 0; i < pso_woa_count; i++)
	    	{
	    		if(i < pso_woa_count/2)
	    		{
	    			woa._whaleSet[i]._fitness = pso_woa_fitness[i];
	    		
	    			for(int j = 0; j < dimension; j++)
	    			{
	    				woa._whaleSet[i]._position[j] = pso_woa_position[pso_woa_index[i]][j];
	    			}
	    		
	    		}
	    		else
	    		{
	    			pso._particleSet[i - pso_woa_count/2]._fitness = pso_woa_fitness[i];
	    		
	    			for(int j = 0; j < dimension; j++)
	    			{
	    				pso._particleSet[i - pso_woa_count/2]._position[j] = pso_woa_position[pso_woa_index[i]][j];
	    			}
	    		}
	    	 
	    	}
	    
	    	woa.refresh();
	    	pso.refresh();
	    	woa.update(Iter, Iter_Max);
	    	pso.update(Iter, Iter_Max, disturbanceRate, disturbanceVelocityCoe);
	    	
	    }
	    
	    for(int i = 0; i < pso_woa_count; i++)
	    {
	    	if(i < pso_woa_count/2)
	    	{
	    		pso_woa_fitness[i] = woa._whaleSet[i]._fitness;
	    		for(int j = 0; j < dimension; j++)
	    		{
	    			pso_woa_position[i][j] = woa._whaleSet[i]._position[j];
	    		}
	    	}
	    	else
	    	{
	    		pso_woa_fitness[i] = pso._particleSet[i - pso_woa_count/2]._fitness;
	    		for(int j = 0; j < dimension; j++)
	    		{
	    			pso_woa_position[i][j] = pso._particleSet[i - pso_woa_count/2]._position[j];
	    		}
	    	}
	    	 pso_woa_index[i] = i;
	    }	            
	    
	    sort_max2min_main(pso_woa_fitness, pso_woa_count, pso_woa_index);
	    
	    cout << pso_woa_fitness[pso_woa_count-1] << endl;
		
		tempt[kk] = pso_woa_fitness[pso_woa_count-1];
		
	}
	    

  
    
    double avgFitness;
    double stddevFitness;
    avgFitness = avg(tempt, num_calc);
    stddevFitness = stddev(tempt, num_calc);
    
    cout << avgFitness << endl;
    cout << stddevFitness << endl;
    

	delete[] minPos;
    delete[] maxPos;
	delete[] tempt;
	delete[] pso_woa_fitness;
	delete[] pso_woa_index;
    delete[] maxSpeed;
    delete[] minSpeed;
	
	for(int i=0;i<pso_woa_count;i++)
		delete[] pso_woa_position[i];

	delete[] pso_woa_position;
    
	
    return(0);
}


//优化目标函数定义
double objFunction_PSO(PSO_Particle& particle)
{    
    double fitnessVal = 0;
    double x_i = 0;
    
    for(int i = 0; i < particle._dimension; i++)
    {
        x_i = particle._position[i];
        
        // fitnessVal += x_i * x_i - (10 * cos(2 * PI * x_i)) + 10;
        // fitnessVal += x_i * x_i;
        fitnessVal += (-(x_i * sin(sqrt(fabs(x_i)))));  // bug
 
    }
        
    return ( fitnessVal );
    
}




double objFunction_WOA(WOA_Whale &whale)
{    
    double fitnessVal = 0;
    double x_i = 0;
    
    for(int i = 0; i < whale._dimension; i++)
    {
        x_i = whale._position[i];
        
        // fitnessVal += x_i * x_i - (10 * cos(2 * PI * x_i)) + 10;
        // fitnessVal += x_i * x_i;
        fitnessVal += (-(x_i * sin(sqrt(fabs(x_i)))));  // bug
        
        
    }
        
    return ( fitnessVal );
    
}


//  4、一些计算指标及随机数函数
// // 求平均数
double avg(double* parameter, int n)
{
	double num = 0;
	for (int i = 0; i < n; ++i)
	{
		num += parameter[i];
	}
	return (num / n);
}
// 求标准差
double stddev(double* parameter, int n)
{
	double num = avg(parameter, n);
	double sum = 0.0;
	for (int i = 0; i < n; ++i)
	{
		sum += (parameter[i] - num) * (parameter[i] - num);
	}
	return sqrt(sum / n);
}


// 数组排序，从大到小排序，并返回排序后的数组对应原数组的下标
void sort_max2min_main(double* a, int length, int* b)
{
    int temp_index;
    double temp;
    for(int j = 0; j < length; j++)
	{
        for(int i = 0; i < length-1-j; i++)
        {
            if(a[i] < a[i+1])
            {
                temp = a[i];
                a[i] = a[i+1];
                a[i+1] = temp;
 
                temp_index = b[i];
                b[i] = b[i+1];
                b[i+1] = temp_index;
            }
        }
	}
}


















