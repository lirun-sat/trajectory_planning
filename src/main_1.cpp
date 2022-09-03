
#include <cmath>
#include <iostream>
#include "TLBO_Algorithm.h"
#include "WOAAlgorithm.h"

using namespace std;


// 采用 经典 函数 测试 woa_tlbo 算法 
//待优化的目标函数
double objFunction_TLBO(Student& student);
double objFunction_WOA(WOA_Whale& whale);
void sort_max2min_main(double* , int , int* );


int main(void)
{         
    double* minPos = new double[30];
    double* maxPos = new double[30];

    for(int i = 0; i < 30; i++)
    {
    	minPos[i] = -500;
    	maxPos[i] = 500;
    }
    
    // int dimension = 3; 
    int dimension = 30;
    int student_whale_count = 40;
    int Iter_Max = 2000;
    int num_calc = 15;

    double* tempt = new double[num_calc];  

	double* woa_tlbo_fitness = new double[student_whale_count];

	int* woa_tlbo_index = new int[student_whale_count];

	double** woa_tlbo_position = new double*[student_whale_count];
	for( int i = 0; i < student_whale_count; i++)
	{
		woa_tlbo_position[i] = new double[dimension];
	}
    
    //构建算法
    WOA_Algorithm woa(objFunction_WOA, minPos, maxPos, dimension, student_whale_count/2);
	
	TLBO_Algorithm tlbo(objFunction_TLBO, minPos, maxPos, dimension, student_whale_count/2);
	
	// *********************************************************************************************************
	// cout << "构建成功" << endl;
	
	for(int kk = 0; kk < num_calc; kk++)
	{
		woa.randomlyInitial();
		tlbo.randomlyInitial();
		
	    for(int iter = 0; iter < Iter_Max; iter++)
	    {
	    	for(int i = 0; i < student_whale_count; i++)
	    	{
	    		if(i < student_whale_count/2)
	    		{
	    			woa_tlbo_fitness[i] = woa._whaleSet[i]._fitness;
	    			for(int j = 0; j < dimension; j++)
	    			{
	    				woa_tlbo_position[i][j] = woa._whaleSet[i]._position[j];
	    			}
	    		}
	    		else
	    		{
	    			woa_tlbo_fitness[i] = tlbo._studentSet[i - student_whale_count/2]._fitness;
	    			for(int j = 0; j < dimension; j++)
	    			{
	    				woa_tlbo_position[i][j] = tlbo._studentSet[i - student_whale_count/2]._position[j];
	    			}
	    		}
	    	}
	    	for(int k = 0; k < student_whale_count; k++)
			{
				woa_tlbo_index[k] = k;
			}
            	 
	    	sort_max2min_main(woa_tlbo_fitness, student_whale_count, woa_tlbo_index);
	    
	    	for(int i = 0; i < student_whale_count; i++)
	    	{
	    		if(i < student_whale_count/2)
	    		{
	    			woa._whaleSet[i]._fitness = woa_tlbo_fitness[i];
	    		
	    			for(int j = 0; j < dimension; j++)
	    			{
	    				woa._whaleSet[i]._position[j] = woa_tlbo_position[woa_tlbo_index[i]][j];
	    			}
	    		}
	    		else
	    		{
	    			tlbo._studentSet[i - student_whale_count/2]._fitness = woa_tlbo_fitness[i];
	    		
	    			for(int j = 0; j < dimension; j++)
	    			{
	    				tlbo._studentSet[i - student_whale_count/2]._position[j] = woa_tlbo_position[woa_tlbo_index[i]][j];
	    			}
	    		}
	    	}
	    
	    	woa.refresh();
	    	tlbo.refresh();
	    	woa.update(iter, Iter_Max);
	    	tlbo.update();
	    	
	    }
	    
	    for(int i = 0; i < student_whale_count; i++)
	    {
	    	if(i < student_whale_count/2)
	    	{
	    		woa_tlbo_fitness[i] = woa._whaleSet[i]._fitness;
	    		for(int j = 0; j < dimension; j++)
	    		{
	    			woa_tlbo_position[i][j] = woa._whaleSet[i]._position[j];
	    		}
	    	}
	    	else
	    	{
	    		woa_tlbo_fitness[i] = tlbo._studentSet[i - student_whale_count/2]._fitness;
	    		for(int j = 0; j < dimension; j++)
	    		{
	    			woa_tlbo_position[i][j] = tlbo._studentSet[i - student_whale_count/2]._position[j];
	    		}
	    	}
	    }	
	    
	    for(int k = 0; k < student_whale_count; k++)
		{
			woa_tlbo_index[k] = k;
		}

	    sort_max2min_main(woa_tlbo_fitness, student_whale_count, woa_tlbo_index);
	    
	    cout << woa_tlbo_fitness[student_whale_count-1] << endl;
		
		tempt[kk] = woa_tlbo_fitness[student_whale_count-1];
		
	}

    double avgFitness;
    double stddevFitness;
    avgFitness = avg(tempt, num_calc);
    stddevFitness = stddev(tempt, num_calc);
    
    cout << "avgFitness: " << avgFitness << endl;
    cout << "stddevFitness: " << stddevFitness << endl;
    
	delete[] minPos;
	delete[] maxPos;
	delete[] tempt;
	delete[] woa_tlbo_fitness;
	delete[] woa_tlbo_index;
	
	for(int i=0;i<student_whale_count;i++)
	{
		delete[] woa_tlbo_position[i];
	}
	delete[] woa_tlbo_position;
    
	
    return(0);
}


//优化目标函数定义
double objFunction_TLBO(Student &student)
{    
    double fitnessVal = 0;
    double x_i = 0;
    
    for(int i = 0; i < student._dimension; i++)
    {
        x_i = student._position[i];
        
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


















