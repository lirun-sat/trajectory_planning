#ifndef _TLBOALGORITHM_H
#define _TLBOALGORITHM_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iomanip>

#define EPS 0.000001


double avg(double* parameter, int n);
double stddev(double* parameter, int n);
void MatrixMultiMatrix(int row_1, int col_1, int col_2, double* A, double* B, double* AB);
int primerange(int start, int end, int* result, int count);
void gen_good_point_set(int num_particles, int dim, double* low_bound, double* up_bound, double* good_point_set);
void bubbleSort(double* arr, int n);


// 算法例子个体
class Student
{
public:
	int _dimension;

	double* _position;

	double _fitness;

	//构造函数，维度初始化为0
	Student(void)
	{
		_dimension = 0;
	}

	//析构函数，释放内存
	~Student(void)
	{
		if (_dimension)
		{
			delete[] _position;
		}
	}

	//初始化函数，用于开辟内存空间
	void initial(int dimension)
	{
		if (_dimension != dimension && dimension)
		{
			//需要重新分配内存
			if (_dimension)
			{
				//消除已有内存
				delete[] _position;
			}
			//开辟新内存
			_dimension = dimension;
			_position = new double[_dimension];
		}
	}
	//复制函数，用于复制操作
	void copy(Student& student)
	{
		this->initial(student._dimension);
		for (int i = 0; i < _dimension; i++)
		{
			_position[i] = student._position[i];
		}
		_fitness = student._fitness;

	}

};


//TLBO算法
class TLBO_Algorithm
{
public:
	int _dimension;
	int _studentCount;
	double _globalBestFitness = 0;
	int _globalBestStudentIndex = -1;

	double r1;
	double r2;
	int TF;  // TF is the teaching factor and is either 1 or 2 (chosen randomly)
	int q;

	Student _best_student;
	Student* _studentSet;

	double *_positionMinValue;
	double *_positionMaxValue;
	double(*_fitnessFunction)(Student&);

	TLBO_Algorithm(double(*objFunction)(Student&), double *positionMinValue, double *positionMaxValue, int dimension, int studentCount)
	{
		//初始化类内参数并分配内存
		_fitnessFunction = objFunction;
		_dimension = dimension;
		_positionMinValue = new double[_dimension];
		_positionMaxValue = new double[_dimension];

		for (int i = 0; i < _dimension; i++)
		{
			_positionMinValue[i] = positionMinValue[i];
			_positionMaxValue[i] = positionMaxValue[i];
		}
		_studentCount = studentCount;
		_studentSet = new Student[_studentCount];

		for (int i = 0; i < _studentCount; i++)
		{
			_studentSet[i].initial(_dimension);
		}
		_best_student.initial(_dimension);

		//配置随机数种子
		srand((unsigned int)time(NULL));

	}
	/***************************************************************
	* 函数描述：析构，释放算法内存
	***************************************************************/
	~TLBO_Algorithm(void)
	{
		//释放内存
		delete[] _positionMinValue;
		delete[] _positionMaxValue;
		delete[] _studentSet;
	}

	double rand0_1(void)
	{
		return((1.0 * rand()) / RAND_MAX);
	}

	// 选出最小的适应度函数值，找到对应的粒子，将其作为最优粒子
	void refresh(void)
	{
		_globalBestStudentIndex = -1;
		double* fitness_temp = new double[_studentCount];

		for (int i = 0; i < _studentCount; i++)
		{
			fitness_temp[i] = _studentSet[i]._fitness;
		}

		bubbleSort(fitness_temp, _studentCount);

		// 选出适应度函数最小的值
		_globalBestFitness = fitness_temp[_studentCount - 1];

		for (int i = 0; i < _studentCount; i++)
		{
			if ((_studentSet[i]._fitness - _globalBestFitness) < EPS)
			{
				_globalBestStudentIndex = i;
			}
			else
			{
				continue;
			}
		}

		if (_globalBestStudentIndex != -1)
		{
			_best_student.copy(_studentSet[_globalBestStudentIndex]);
		}
			
		delete[] fitness_temp;

	}


	void randomlyInitial(void)
	{
		double* good_point_set_temp = new double[_studentCount* _dimension];

		//  产生佳点集
		gen_good_point_set(_studentCount, _dimension, _positionMinValue, _positionMaxValue, good_point_set_temp);
		
		for (int i = 0; i < _studentCount; i++)
		{
			for (int j = 0; j < _dimension; j++)
			{
				_studentSet[i]._position[j] = good_point_set_temp[i * _dimension + j];
			}
			_studentSet[i]._fitness = _fitnessFunction(_studentSet[i]);
		}

		delete[] good_point_set_temp;

	}


	void update(void)
	{
		double* X_mean = new double[_dimension];
		double X_mean_temp;

		Student Xnew_temp_1;
		Xnew_temp_1.initial(_dimension);

		Student Xnew_temp_2;
		Xnew_temp_2.initial(_dimension);

		for (int i = 0; i < _studentCount; i++)
		{
			for (int j = 0; j < _dimension; j++)
			{
				X_mean_temp = 0;
				for (int k = 0; k < _studentCount; k++)
				{
					X_mean_temp += _studentSet[k]._position[j];
				}
				X_mean[j] = X_mean_temp / _studentCount;
			}

			TF = (rand() % (2 - 1 + 1)) + 1;   //  要取得[a,b]的随机整数，使用(rand() % (b-a+1))+ a;

			// Teaching phase-----------
			for (int j = 0; j < _dimension; j++)
			{
				r1 = rand0_1();

				Xnew_temp_1._position[j] = _studentSet[i]._position[j] + r1 * (_best_student._position[j] - TF * X_mean[j]);

				if (Xnew_temp_1._position[j] >= _positionMaxValue[j])
				{
					Xnew_temp_1._position[j] = _positionMaxValue[j];
				}
				else if (Xnew_temp_1._position[j] <= _positionMinValue[j])
				{
					Xnew_temp_1._position[j] = _positionMinValue[j];
				}
			}

			Xnew_temp_1._fitness = this->_fitnessFunction(Xnew_temp_1);

			if (Xnew_temp_1._fitness < _studentSet[i]._fitness)
			{
				for (int j = 0; j < _dimension; j++)
				{
					_studentSet[i]._position[j] = Xnew_temp_1._position[j];
				}
				_studentSet[i]._fitness = Xnew_temp_1._fitness;
			}

			if (Xnew_temp_1._fitness < _best_student._fitness)
			{
				for (int j = 0; j < _dimension; j++)
				{
					_best_student._position[j] = Xnew_temp_1._position[j];
				}
				_best_student._fitness = Xnew_temp_1._fitness;
			}


			// Learning phase----------------------------------------------------------------------------------------------------------------------------------------------
			q = (rand() % ((_studentCount - 1) - 0 + 1)) + 0;   //  要取得[a,b]的随机整数，使用(rand() % (b-a+1))+ a;

			while (q == i)
			{
				q = (rand() % ((_studentCount - 1) - 0 + 1)) + 0;
			}

			if (_studentSet[i]._fitness < _studentSet[q]._fitness)
			{
				for (int j = 0; j < _dimension; j++)
				{
					r2 = rand0_1();
					Xnew_temp_2._position[j] = _studentSet[i]._position[j] + r2 * (_studentSet[i]._position[j] - _studentSet[q]._position[j]);
				}
			}
			else
			{
				for (int j = 0; j < _dimension; j++)
				{
					r2 = rand0_1();
					Xnew_temp_2._position[j] = _studentSet[i]._position[j] - r2 * (_studentSet[i]._position[j] - _studentSet[q]._position[j]);
				}
			}

			for (int j = 0; j < _dimension; j++)
			{
				if (Xnew_temp_2._position[j] >= _positionMaxValue[j])
				{
					Xnew_temp_2._position[j] = _positionMaxValue[j];
				}
				else if (Xnew_temp_2._position[j] <= _positionMinValue[j])
				{
					Xnew_temp_2._position[j] = _positionMinValue[j];
				}
			}

			Xnew_temp_2._fitness = this->_fitnessFunction(Xnew_temp_2);

			if (Xnew_temp_2._fitness < _studentSet[i]._fitness)
			{
				for (int j = 0; j < _dimension; j++)
				{
					_studentSet[i]._position[j] = Xnew_temp_2._position[j];
				}
				_studentSet[i]._fitness = Xnew_temp_2._fitness;
			}

			if (Xnew_temp_2._fitness < _best_student._fitness)
			{
				for (int j = 0; j < _dimension; j++)
				{
					_best_student._position[j] = Xnew_temp_2._position[j];
				}
				_best_student._fitness = Xnew_temp_2._fitness;
			}
		}

		delete[] X_mean;

	}


	// void findMin(int max_iter, Student& bestStudent)
	// {
	// 	this->randomlyInitial();

	// 	for (int iter = 0; iter < max_iter; iter++)
	// 	{
	// 		this->update();
	// 	}

	// 	bestStudent.copy(_best_student);
	// }

};


// 仅仅找出最小的元素，并将其放在最后
void bubbleSort(double* arr, int n)
{
	double temp;
	//比较两个相邻的元素   
	for (int j = 0; j < n - 1; j++)
	{
		if (arr[j] < arr[j + 1])
		{
			temp = arr[j];
			arr[j] = arr[j + 1];
			arr[j + 1] = temp;
		}
	}
}

//  4、一些计算指标及随机数函数
double avg(double* parameter, int n)
{
    // 求平均数
	double num = 0;
	for (int i = 0; i < n; i++) 
	{
		num += parameter[i];
	}
    
	return (num / n);
    
}


double stddev(double* parameter, int n) 
{
    // 求标准差
	double num = avg(parameter, n);
	double sum = 0.0;
    
	for (int i = 0; i < n; i++) 
	{
		sum += (parameter[i] - num) * (parameter[i] - num);
	}
    
	return sqrt(sum / n);
    
}



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





#endif 

