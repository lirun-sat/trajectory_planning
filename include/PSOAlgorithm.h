#ifndef _ZPSOALGORITHM_H
#define _ZPSOALGORITHM_H

#define _USE_MATH_DEFINES 

#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <random>
#include <chrono>



#define WMAX 0.9                           
#define WMIN 0.1      

void MatrixMultiMatrix(int row_1, int col_1, int col_2, double* A, double* B, double* AB);
int primerange(int start, int end, int* result, int count);
void gen_good_point_set(int num_particles, int dim, double* low_bound, double* up_bound, double* good_point_set);
void bubbleSort(double* arr, int n);


class PSO_Particle
{
public:
    int _dimension;        
    double *_position;     
    double *_velocity;     
    double *_bestPosition; 
    double _fitness;       
    double _bestFitness;   
    
    PSO_Particle(void){
        _dimension = 0;
    }
    
    ~PSO_Particle(void){
        if(_dimension){
            delete [] _position;
            delete [] _velocity;
            delete [] _bestPosition;
        }
    }
    
    void initial(int dimension){
        if(_dimension!=dimension && dimension){
            if(_dimension){
                //消除已有内存
                delete [] _position;
                delete [] _velocity;
                delete [] _bestPosition;
            }
            //开辟新内存
            _dimension = dimension;
            _position = new double[_dimension];
            _velocity = new double[_dimension];
            _bestPosition = new double[_dimension];
            
        }
        
    }
    //复制函数，用于粒子间的复制操作
    void copy(PSO_Particle& particle){
        this->initial(particle._dimension);
        for(int i = 0; i < _dimension; i++){
            _position[i] = particle._position[i];
            _velocity[i] = particle._velocity[i];
            _bestPosition[i] = particle._bestPosition[i];
            _fitness = particle._fitness;
            _bestFitness = particle._bestFitness;
            
        }
    }
};


class PSO_Algorithm
{
public:
    int _dimension;         // 粒子群维度
    int _particleCount;     // 种群粒子数量
    double _globalGuideCoe; // 全局最优引导系数
    double _localGuideCoe;  // 局部最优引导系数
    double _inertGuideCoe;  // inertia coefficient 
    double _globalBestParticleFitness = 0;  
    double* _positionMinValue;  // 粒子位置的最小界
    double* _positionMaxValue;  // 粒子位置的最大界
    double* _maxSpeed;          // 粒子允许最大速度
    double* _minSpeed;
    double (*_fitnessFunction)(PSO_Particle&);  //粒子适应度函数
    PSO_Particle* _particleSet;  // 粒子集
    PSO_Particle _globalBestParticle;    // 搜索过程得到的全局最优粒子

    PSO_Algorithm(double (*objFunction)(PSO_Particle&), double* positionMinValue, double* positionMaxValue, 
                  int dimension, int particleCount, double* maxSpeed, double* minSpeed, 
                  double inertGuideCoe, double globalGuideCoe = 2, double localGuideCoe = 2){
        //初始化类内参数并分配内存
        _fitnessFunction = objFunction;
        _dimension = dimension;
        _positionMinValue = new double[_dimension];
        _positionMaxValue = new double[_dimension];
        _maxSpeed = new double[_dimension];
        _minSpeed = new double[_dimension];
        _particleSet = new PSO_Particle[_particleCount];

        for(int i = 0; i < _dimension; i++){
            _positionMinValue[i] = positionMinValue[i];
            _positionMaxValue[i] = positionMaxValue[i];
            _maxSpeed[i] = maxSpeed[i];
            _minSpeed[i] = minSpeed[i];
        }
        _particleCount = particleCount;
        _globalGuideCoe = globalGuideCoe;
        _localGuideCoe = localGuideCoe;
        _inertGuideCoe = inertGuideCoe;

        for(int i = 0; i< _particleCount; i++){
            _particleSet[i].initial(_dimension);
        }   
        _globalBestParticle.initial(_dimension);
        
        //配置随机数种子
        srand((unsigned int)time(NULL));
        
    }
    /***************************************************************
     * 函数描述：析构一个PSO算法，释放算法内存
    ***************************************************************/
    ~PSO_Algorithm(void){
        //释放内存
        delete [] _positionMinValue;
        delete [] _positionMaxValue;
        delete [] _particleSet;
        delete [] _maxSpeed;
        delete [] _minSpeed;
    }

    double rand0_1(void){
        return((1.0 * rand()) / RAND_MAX);
    }

    void randomlyInitial(void){
        int globalBestParticleIndex = -1;
        double* good_point_set_temp = new double[_particleCount * _dimension];
        double* fitness_temp = new double[_particleCount];
        
        //  产生佳点集
        gen_good_point_set(_particleCount, _dimension, _positionMinValue, _positionMaxValue, good_point_set_temp);

        for(int i = 0; i < _particleCount; i++){
            for(int j = 0; j < _dimension; j++){
                _particleSet[i]._position[j] = good_point_set_temp[i * _dimension + j];
                _particleSet[i]._bestPosition[j] = _particleSet[i]._position[j];
                if(rand0_1() > 0.5){
                    _particleSet[i]._velocity[j] = rand0_1() * _maxSpeed[j];
                }
                else{
                    _particleSet[i]._velocity[j] = rand0_1() * _minSpeed[j];
                }
            }
            _particleSet[i]._fitness = _fitnessFunction(_particleSet[i]);
            _particleSet[i]._bestFitness = _particleSet[i]._fitness;
            fitness_temp[i] = _particleSet[i]._bestFitness;
        }
        bubbleSort(fitness_temp, _particleCount);
        _globalBestParticleFitness = fitness_temp[ _particleCount - 1];
        for(int i = 0; i < _particleCount; i++){
            if(_particleSet[i]._bestFitness <= _globalBestParticleFitness){
                globalBestParticleIndex = i;
            }
        }
        if(globalBestParticleIndex != -1){
            _globalBestParticle.copy(_particleSet[globalBestParticleIndex]);
        }
            
        delete[] good_point_set_temp;
        delete[] fitness_temp;
       
    }

    void refresh(void){
        int globalBestParticleIndex = -1;
        for(int i = 0; i < _particleCount; i++){
            _particleSet[i]._fitness = this->_fitnessFunction(_particleSet[i]);  
            if(_particleSet[i]._fitness < _particleSet[i]._bestFitness){
                for(int j = 0; j < _dimension; j++){
                    _particleSet[i]._bestPosition[j] = _particleSet[i]._position[j];
                }
                _particleSet[i]._bestFitness = _particleSet[i]._fitness;
                //是否更新全局最优解
                if(_particleSet[i]._bestFitness < _globalBestParticleFitness){
                    globalBestParticleIndex = i;
                    _globalBestParticleFitness = _particleSet[i]._bestFitness;  
                }
            }
        }
        //更新全局最优粒子位置
        if(globalBestParticleIndex != -1){
            _globalBestParticle.copy(_particleSet[globalBestParticleIndex]);
        }
    }
    
    /***************************************************************
     * 函数名：disturbance
     * 函数描述：对粒子速度进行给定大小的扰动
     * 输入参数：
     *  particle：被扰动的粒子对象
     *  relativeVelocityRate：扰动速度大小上限相对于_maxSpeed的比例，默认为0.05
     * 输出参数：void
    ***************************************************************/
    void disturbance(PSO_Particle &particle, double disturbanceVelocityCoe){
        // 生成扰动速度
        double* disturbanceVelocity = new double[_dimension];
        
        for(int i = 0; i < _dimension; i++){
            if(rand0_1() > 0.5){
                disturbanceVelocity[i] = rand0_1() * disturbanceVelocityCoe * _maxSpeed[i];
            }
            else{
                disturbanceVelocity[i] = rand0_1() * disturbanceVelocityCoe * _minSpeed[i];
            }

            particle._velocity[i] += disturbanceVelocity[i];
            
            if(particle._velocity[i] > _maxSpeed[i]){
                particle._velocity[i] = _maxSpeed[i];
            }
            else if(particle._velocity[i] < _minSpeed[i]){
                particle._velocity[i] = _minSpeed[i];
            }
        }
                
		delete[] disturbanceVelocity;
		
    }

    void update(int iter, int Iter_Max, double disturbanceRate, double disturbanceVelocityCoe){
        _inertGuideCoe = WMAX - (WMAX - WMIN) / Iter_Max * (double)iter;  //  需要对 i 进行类型转换， 因为 _inertGuideCoe 是 double 型变量.

        for(int i = 0; i < _particleCount; i++){
            double r1 = rand0_1();
            double r2 = rand0_1();
            for(int j = 0; j < _particleSet[i]._dimension; j++){
                _particleSet[i]._velocity[j] = _inertGuideCoe * _particleSet[i]._velocity[j] 
                                                + _globalGuideCoe * r1 * (_globalBestParticle._bestPosition[j] - _particleSet[i]._position[j]) 
                                                + _localGuideCoe * r2 * (_particleSet[i]._bestPosition[j] - _particleSet[i]._position[j]); 

                if(_particleSet[i]._velocity[j] > _maxSpeed[j]){
                    _particleSet[i]._velocity[j] = _maxSpeed[j];
                } 
                else if(_particleSet[i]._velocity[j] < _minSpeed[j]){
                    _particleSet[i]._velocity[j] = _minSpeed[j];
                }                               
            }
            //对粒子速度进行扰动，提高算法局部搜索能力
            if(rand0_1() < disturbanceRate){
                this->disturbance(_particleSet[i], disturbanceVelocityCoe);
            }
            for(int j = 0; j < _particleSet[i]._dimension; j++){
                _particleSet[i]._position[j] += _particleSet[i]._velocity[j];
                if(_particleSet[i]._position[j] < _positionMinValue[j]){
                    _particleSet[i]._position[j] = _positionMinValue[j];
                }
                else if(_particleSet[i]._position[j] > _positionMaxValue[j]){
                    _particleSet[i]._position[j] = _positionMaxValue[j];
                }
            }
        }
    }
    
};


void MatrixMultiMatrix(int row_1, int col_1, int col_2, double* A, double* B, double* AB)
{
	
	/*
	 两个 任意维数 矩阵相乘
	 row_1 第一个矩阵的行数
	 col_1 第一个矩阵的列数
	 col_1 第二个矩阵的行数和第一个矩阵的列数相等
	 col_2 第二个矩阵的列数
	*/
    
    double* A_tempt;
    A_tempt = new double[row_1 * col_1];
    
    double* B_tempt;
    B_tempt = new double[col_1 * col_2];
    
    double* AB_tempt;
    AB_tempt = new double[row_1 * col_2];
    
    for(int i = 0; i < row_1; i++)
    {
    	for(int j = 0; j < col_1; j++)
    	{
    		A_tempt[i * col_1 + j] = A[i * col_1 + j];
    	}
    }
    
    for(int i = 0; i < col_1; i++)
    {	
        for(int j = 0; j < col_2; j++)
        {
    		B_tempt[i * col_2 + j] = B[i * col_2 + j];
    	}
    }
    
    for(int i = 0; i < row_1; i++)
    {
        for(int j = 0; j < col_2; j++)
        {
    	    AB_tempt[i * col_2 + j] = 0;
    	}
    }
    		
	for(int i = 0; i < row_1; i++)
	{
		for(int j = 0; j < col_2; j++)
		{
			for(int k = 0; k < col_1; k++)
			{
				AB_tempt[i * col_2 + j] += A_tempt[i * col_1 + k] * B_tempt[k * col_2 + j];
			}
		}
	}
    
    for(int i = 0; i < row_1; i++)
    {
    	for(int j = 0; j < col_2; j++)
    	{
    		AB[i * col_2 + j] = AB_tempt[i * col_2 + j];
    	}
    }
    
    delete[] A_tempt;
    delete[] B_tempt;
    delete[] AB_tempt;
   	
}


int primerange(int start, int end, int* result, int count)
{
    // 求 [start, end] 之间的 素数， 返回值 count 表示 素数总个数 
    int* result_temp = new int[1000];
    int result_count = 0;
    bool prime;  //定义bool变量
            
    
    for(int i = start; i < (end+1); i++)
    {
        prime = true;    //先令 prime 为真
        
        if((i == 0) || (i == 1))
        {
            // cout << "i=0 or i=1" <<endl;
            continue;
        }
       
        for(int j = 2; j < i; j++)    //对 2 到 m 进行循环
        {
            if(i % j == 0)    //若 i 整除 j 为 0，令 prime 为假，循环终止
            {
                prime = false;
                break;
            }
        }
        if(prime)    //若 prime 为真，输出 n
        {
            result_temp[result_count] = i;
            result_count += 1;
        }
        
    }
    
    for(int i = 0; i < result_count; i++)
        result[i] = result_temp[i];
    
    count = result_count;
    
    delete[] result_temp;
    
    return count;
    
}


void bubbleSort(double* arr, int n)
{
    //  进行排序, max -> min
    double temp;
    
	for(int i = 0; i < n; i++)
	{  
		//比较两个相邻的元素   
		for(int j = 0; j < n-i-1; j++)
		{  
            if(arr[j] < arr[j+1])
            {  
                temp = arr[j];  
                arr[j] = arr[j+1];  
                arr[j+1] = temp;  
            }  
        }  
    }       
}



void gen_good_point_set(int num_particles, int dim, double* low_bound, double* up_bound, double* good_point_set)
{
    // num_particles: 群体个数； dim: 每个个体的维数； low_bound: 各个维度的 下限值；up_bound: 各个维度的 上限值； good_point_set: 生成的 佳点集
    double pi = 3.14159265358979323846;
    
    double* temp_1 = new double[num_particles];
    double* ones_dim = new double[dim];
    double* temp_1_multi_ones_dim = new double[num_particles * dim];
    double* up_bound_minus_low_bound = new double[dim];
    double* ind = new double[dim];
    double* temp_2 = new double[dim];
    double* ones_num_particles = new double[num_particles];
    double* ones_num_particles_multi_temp_2 = new double[num_particles * dim];
    double* good_points = new double[num_particles * dim];
    double* good_points_temp = new double[num_particles * dim];
    int prime_count = 0;
    int* prime_list_tempt = new int[1000];
    // int idx;
    int prime_idx = 0;
    
    
    for(int i=0;i<num_particles;i++)
    {
        temp_1[i] = i+1;
        
    }
    for(int i=0;i<dim;i++)
    {
        ones_dim[i] = 1;
    }
    
    MatrixMultiMatrix( num_particles,  1,  dim, temp_1, ones_dim, temp_1_multi_ones_dim);
    
    // ind = np.arange(1, dim + 1)
    for(int i=0;i<dim;i++)
        ind[i] = i+1;
    
    // prime = list(sympy.sieve.primerange(0, 100 * dim))
    prime_count = primerange( 0, 100 * dim, prime_list_tempt, prime_count);
    
    int* prime_list = new int[prime_count];    //  注意，这是定义 prime_list，不能提前定义！
    for(int i=0;i<prime_count;i++)
        prime_list[i] = prime_list_tempt[i];
    
    for(int i=0;i<prime_count;i++)
    {
        if(prime_list[i] < (2 * dim + 3) )
            continue;
        
        else
        {
            // idx = i;
            prime_idx = prime_list[i];
            break;
        }
        
    }
    
    for(int i=0;i<dim;i++)
    {
        temp_2[i] = 2 * pi * ind[i] / prime_idx;
       
        temp_2[i] = 2 * cos(temp_2[i]);
        
    }
    

    for(int i = 0; i < num_particles; i++)
    {
        ones_num_particles[i] = 1;
    }
    
    MatrixMultiMatrix( num_particles, 1, dim, ones_num_particles, temp_2, ones_num_particles_multi_temp_2);
    
    //  gd = temp1 * temp2
    for(int i = 0; i < (num_particles * dim); i++)
        good_points[i] = temp_1_multi_ones_dim[i] * ones_num_particles_multi_temp_2[i];
    
    for(int i = 0; i < (num_particles * dim); i++)
        good_points_temp[i] = good_points[i];
    
    // gd = np.mod(gd, 1)
    for(int i = 0; i < (num_particles * dim); i++)
    {
        good_points_temp[i] = fmod(good_points_temp[i], 1);
        
        if(good_points_temp[i] < 0)
            good_points_temp[i] = 1 + good_points_temp[i];
    
    }
    
    for(int i = 0; i < (num_particles * dim); i++)
        good_points[i] = good_points_temp[i];
    
    /*
    aa = list(map(lambda x, y: x - y, up_bound, low_bound))
    */
    for(int i = 0; i < dim; i++)
    {
        up_bound_minus_low_bound[i] = up_bound[i] - low_bound[i];
    }
    
    
    for(int i = 0; i < num_particles; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            //  bb = list(map(lambda x, y: x * y, [gd[:, i] for i in range(dim)], aa))
            good_points[i*dim + j] = good_points[i*dim + j] * up_bound_minus_low_bound[j];
            
            //  cc = list(map(lambda x, y: x + y, bb, low_bound))
            good_points[i*dim + j] = good_points[i*dim + j] + low_bound[j];
        }
    }
    
    for(int i = 0; i < (num_particles * dim); i++)
        good_point_set[i] = good_points[i];
    
    
    delete[] temp_1;
    delete[] ones_dim;
    delete[] temp_1_multi_ones_dim;
    delete[] up_bound_minus_low_bound;
    delete[] ind;
    delete[] temp_2;
    delete[] ones_num_particles;
    delete[] ones_num_particles_multi_temp_2;
    delete[] good_points;
    delete[] good_points_temp;
    delete[] prime_list_tempt;
    delete[] prime_list;
     
}




#endif // _ZPSOALGORITHM_H


