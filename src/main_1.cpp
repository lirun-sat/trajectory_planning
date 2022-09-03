#include "TLBO_Algorithm.h"
#include <cmath>
#include <iostream>

using namespace std;


// 待优化的目标函数
double objFunction(Student& student);


int main(void)
{
    // 定义求解范围,x在[-10,10]范围，y在[-10,10]范围
    // double minPos[] = {-10, -10, -10};
    // double maxPos[] = {10, 10, 10};
    
    double minPos[] = {-500, -500, -500};
    double maxPos[] = {500, 500, 500};
    
    // 定义问题描述参数
    int dimension = 3; 
    int studentCount = 50; 
    int iteration_max = 800; 
    int num_calc = 20;
    double* tempt = new double[num_calc];  //  计算 20 次， 求 计算的 均值 和 标准差
    
    // 构建算法
    TLBO_Algorithm tlbo(objFunction, minPos, maxPos, dimension, studentCount);

    Student bestStudent;    // 最终进化结果
    bestStudent.initial(dimension); 
    
    for(int i = 0; i < num_calc; i++)
    {
    	// 获取最终进化结果
        tlbo.findMin(iteration_max, bestStudent);    
        tempt[i] = bestStudent._fitness;
        cout <<  tempt[i] << endl;
    }
    
    double avgFitness;
    double stddevFitness;
    avgFitness = avg(tempt, num_calc);
    stddevFitness = stddev(tempt, num_calc);
    
    cout << "avgFitness:" << avgFitness << endl;
    cout << "stddevFitness:" << stddevFitness << endl;
    
    delete[] tempt;
    
    return(0);
}


//优化目标函数定义
double objFunction(Student &student)
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




