#include "WOAAlgorithm.h"
#include <cmath>
#include <iostream>

using namespace std;


//待优化的目标函数
double objFunction(WOA_Whale& whale);


int main(void)
{
    //定义求解范围,x在[-10,10]范围，y在[-10,10]范围
    // double minPos[] = {-500, -500, -500};
    // double maxPos[] = {500, 500, 500};
    
    // //定义问题描述参数
    // int dimension = 3; 

    double* minPos = new double[30];
    double* maxPos = new double[30];

    for(int i = 0; i < 30; i++)
    {
    	minPos[i] = -500;
    	maxPos[i] = 500;
    }
    
    // int dimension = 3; 
    int dimension = 30;
    int whaleCount = 30; 
    int iteration_max = 3000;
    int num_calc = 20;
    double* tempt = new double[num_calc];  //  计算 20 次， 求 计算的 均值 和 标准差
    
    //构建gwo算法
    WOA_Algorithm woa(objFunction, minPos, maxPos, dimension, whaleCount);

    WOA_Whale bestWhale;  //最终进化结果
    
    for(int i = 0; i < num_calc; i++)
    {
    	// TLBO_Algorithm tlbo(objFunction, minPos, maxPos, dimension, studentCount);
    	// srand((unsigned int)time(NULL));
        woa.findMin(iteration_max, bestWhale);  //获取最终进化结果
        tempt[i] = bestWhale._fitness;
        cout <<  tempt[i] << endl;
    }
    
    double avgFitness;
    double stddevFitness;
    avgFitness = avg(tempt, num_calc);
    stddevFitness = stddev(tempt, num_calc);
    
    cout << "avgFitness: " << avgFitness << endl;
    cout << "stddevFitness: " << stddevFitness << endl;
    
    
    delete[] tempt;
    delete[] minPos;
    delete[] maxPos;
    
    return(0);
}


//优化目标函数定义
double objFunction(WOA_Whale &whale)
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
        
    return( fitnessVal );
    
}



