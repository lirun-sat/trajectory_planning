#include "WOAAlgorithm.h"
#include <cmath>
#include <iostream>

using namespace std;


//待优化的目标函数
double objFunction(WOA_Whale& whale);


int main(void)
{
    //定义求解范围,x在[-10,10]范围，y在[-10,10]范围
    double minPos[] = {-500, -500, -500};
    double maxPos[] = {500, 500, 500};
    
    //定义问题描述参数
    int dimension = 3; 
    int whaleCount = 30; 
    int generation = 800; //粒子群进化代数
    int num_calc = 10;
    double* tempt;
    tempt = new double[num_calc];  //  计算 200 次， 求 计算的 均值 和 标准差
    
    //构建gwo算法
    WOA_Algorithm woa(objFunction, minPos, maxPos, dimension, whaleCount);

    WOA_Whale bestWhale; //粒子群最终进化结果
    
    for(int i = 0; i < num_calc; i++)
    {
    	// TLBO_Algorithm tlbo(objFunction, minPos, maxPos, dimension, studentCount);
    	// srand((unsigned int)time(NULL));
        woa.findMin(generation, bestWhale); //获取最终进化结果
        tempt[i] = bestWhale._fitness;
        cout <<  tempt[i] << endl;
    }
    
    double avgFitness;
    double stddevFitness;
    avgFitness = avg(tempt, num_calc);
    stddevFitness = stddev(tempt, num_calc);
    
    cout << avgFitness << endl;
    cout << stddevFitness << endl;
    
    
    // delete[] bestParticleset;
    
    delete[] tempt;
    
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




