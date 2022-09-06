
#ifndef _WOAALGORITHM_H
#define _WOAALGORITHM_H

#define _USE_MATH_DEFINES 

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

void sort_max2min(double* a, int length, int* b);


class WOA_Whale
{
public:
	int _dimension; //粒子维度
	double *_position; //粒子所在位置数组指针
	double _fitness; //例子适应度

					 //构造函数，粒子维度初始化为0
	WOA_Whale(void)
	{
		_dimension = 0;
	}

	//析构函数，释放粒子内存
	~WOA_Whale(void)
	{
		if (_dimension)
		{
			delete[] _position;
		}
	}

	//初始化函数，用于为粒子开辟内存空间
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
	//复制函数，用于粒子间的复制操作
	void copy(WOA_Whale& whale)
	{
		this->initial(whale._dimension);
		for (int i = 0; i<_dimension; i++)
		{
			_position[i] = whale._position[i];
			_fitness = whale._fitness;
		}

	}

};

//PSO算法
class WOA_Algorithm
{
public:
	int _dimension;
	int _whaleCount;

	double a, a2;
	double A, C, b, l;

	double _globalBestFitness = 0;
	int _globalBestWhaleIndex = -1;

	WOA_Whale _best_Whale;

	WOA_Whale* _whaleSet;
	double *_positionMinValue;
	double *_positionMaxValue;
	double(*_fitnessFunction)(WOA_Whale&);

	WOA_Algorithm(double(*objFunction)(WOA_Whale&), double *positionMinValue, double *positionMaxValue, int dimension, int whaleCount)
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

		_whaleCount = whaleCount;
		_whaleSet = new WOA_Whale[_whaleCount];

		for (int i = 0; i < _whaleCount; i++)
			_whaleSet[i].initial(_dimension);

		_best_Whale.initial(_dimension);

		//配置随机数种子
		srand((unsigned int)time(NULL));

	}
	/***************************************************************
	* 函数描述：析构一个PSO算法，释放算法内存
	***************************************************************/
	~WOA_Algorithm(void)
	{
		//释放内存
		delete[] _positionMinValue;
		delete[] _positionMaxValue;
		delete[] _whaleSet;

	}

	double rand0_1(void)
	{
		return((1.0 * rand()) / RAND_MAX);
	}


	// 对 _whaleSet 按照 _fitness 从大到小排列, 并对 _best_Whale 赋值
	void refresh(void)
	{
		_globalBestWhaleIndex = -1;

		double* fitness_temp = new double[_whaleCount];

		WOA_Whale* _whaleSet_temp = new WOA_Whale[_whaleCount];

		int* index = new int[_whaleCount];

		for (int i = 0; i < _whaleCount; i++)
			fitness_temp[i] = _whaleSet[i]._fitness;

		for (int k = 0; k < _whaleCount; k++)
			index[k] = k;

		sort_max2min(fitness_temp, _whaleCount, index);  // fitness_temp 从大到小排序

														 // 选出适应度函数最小的值
		_globalBestFitness = fitness_temp[_whaleCount - 1];

		_globalBestWhaleIndex = index[_whaleCount - 1];

		if (_globalBestWhaleIndex != -1)
		{
			_best_Whale.copy(_whaleSet[_globalBestWhaleIndex]);
		}


		for (int i = 0; i < _whaleCount; i++)
		{
			_whaleSet_temp[i].copy(_whaleSet[index[i]]);
		}

		for (int i = 0; i < _whaleCount; i++)
		{
			_whaleSet[i].copy(_whaleSet_temp[i]);
		}

		delete[] fitness_temp;
		delete[] index;
		delete[] _whaleSet_temp;

	}


	void randomlyInitial(void)
	{
		double* good_point_set_temp;
		good_point_set_temp = new double[_whaleCount * _dimension];
		//  产生佳点集
		gen_good_point_set(_whaleCount, _dimension, _positionMinValue, _positionMaxValue, good_point_set_temp);

		for (int i = 0; i < _whaleCount; i++)
		{
			for (int j = 0; j < _dimension; j++)
				_whaleSet[i]._position[j] = good_point_set_temp[i * _dimension + j];

			_whaleSet[i]._fitness = _fitnessFunction(_whaleSet[i]);
		}

		// this->refresh(); 

		delete[] good_point_set_temp;

	}

	void update(void)
	{
		double* X_D;
		X_D = new double[_dimension];

		double* Xnew;
		Xnew = new double[_dimension];

		double* X_D_prime;
		X_D_prime = new double[_dimension];

		double* Xrand;
		Xrand = new double[_dimension];

		double r1;
		double r2;
		double r3;

		double p;
		int q;

		for (int i = 0; i < _whaleCount; i++)
		{
			r1 = rand0_1();
			r2 = rand0_1();
			r3 = rand0_1();

			b = 1;
			A = a * (2 * r1 - 1);
			C = 2 * r2;
			l = (a2 - 1) * r3 + 1;

			p = rand0_1();

			if (p < 0.5)
			{
				if (fabs(A) > 1)
				{
					q = (rand() % ((_whaleCount - 1) - 0 + 1)) + 0;   //  要取得[a,b]的随机整数，使用(rand() % (b-a+1))+ a;
																	  // q = (rand() % (_whaleCount - 0)) + 0;  // 要取得[a,b)的随机整数，使用(rand() % (b-a))+ a;

					while (q == i)
						q = (rand() % ((_whaleCount - 1) - 0 + 1)) + 0;

					for (int j = 0; j < _dimension; j++)
					{
						Xrand[j] = _whaleSet[q]._position[j];

						X_D[j] = fabs(C * Xrand[j] - _whaleSet[i]._position[j]);
						Xnew[j] = Xrand[j] - A * X_D[j];

					}

				}
				else
				{
					for (int j = 0; j < _dimension; j++)
					{
						X_D[j] = fabs(C * _best_Whale._position[j] - _whaleSet[i]._position[j]);
						Xnew[j] = _best_Whale._position[j] - A * X_D[j];
					}
				}

			}
			else
			{
				for (int j = 0; j < _dimension; j++)
				{
					// X_D_prime[j] = fabs(C * _best_Whale._position[j] - _whaleSet[i]._position[j]);  //  BUG

					X_D_prime[j] = fabs(_best_Whale._position[j] - _whaleSet[i]._position[j]);
					Xnew[j] = X_D_prime[j] * exp(b * l) * cos(2 * M_PI * l) + _best_Whale._position[j];
				}

			}

			for (int j = 0; j < _dimension; j++)
			{
				if (Xnew[j] <= _positionMaxValue[j] && Xnew[j] >= _positionMinValue[j])
					_whaleSet[i]._position[j] = Xnew[j];

				else if (Xnew[j] > _positionMaxValue[j])
					_whaleSet[i]._position[j] = _positionMaxValue[j];

				else
					_whaleSet[i]._position[j] = _positionMinValue[j];

			}

			_whaleSet[i]._fitness = this->_fitnessFunction(_whaleSet[i]);

		}


		// this->refresh();

		delete[] X_D;
		delete[] Xnew;
		delete[] X_D_prime;
		delete[] Xrand;

	}


	void findMin(int max_iter, WOA_Whale& bestWhale)
	{
		this->randomlyInitial();
		this->refresh();

		for (int iter = 0; iter < max_iter; iter++)
		{
			a = 2 * (1 - (double)iter / (double)max_iter);  //  需要对 iter 进行类型转换， 因为 a 是 double 型变量.
			a2 = (-1) + (double)iter * ((-1) / (double)max_iter);

			this->update();
			this->refresh();

		}

		bestWhale.copy(_best_Whale);
	}

};



//  4、一些计算指标及随机数函数
// // 求平均数
double avg(double* parameter, int n)
{
	double num = 0;
	for (int i = 0; i < n; i++)
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

	for (int i = 0; i < n; i++)
	{
		sum += (parameter[i] - num) * (parameter[i] - num);
	}

	return sqrt(sum / n);

}



// 数组排序，从大到小排序，并返回排序后的数组对应原数组的下标
void sort_max2min(double* a, int length, int* b)
{
	int temp_index;
	double temp;
	for (int j = 0; j < length; j++)
	{
		for (int i = 0; i < length - 1 - j; i++)
		{
			if (a[i] < a[i + 1])
			{
				temp = a[i];
				a[i] = a[i + 1];
				a[i + 1] = temp;

				temp_index = b[i];
				b[i] = b[i + 1];
				b[i + 1] = temp_index;
			}
		}
	}
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

	for (int i = 0; i < row_1; i++)
	{
		for (int j = 0; j < col_1; j++)
		{
			A_tempt[i * col_1 + j] = A[i * col_1 + j];
		}
	}

	for (int i = 0; i < col_1; i++)
	{
		for (int j = 0; j < col_2; j++)
		{
			B_tempt[i * col_2 + j] = B[i * col_2 + j];
		}
	}

	for (int i = 0; i < row_1; i++)
	{
		for (int j = 0; j < col_2; j++)
		{
			AB_tempt[i * col_2 + j] = 0;
		}
	}

	for (int i = 0; i < row_1; i++)
	{
		for (int j = 0; j < col_2; j++)
		{
			for (int k = 0; k < col_1; k++)
			{
				AB_tempt[i * col_2 + j] += A_tempt[i * col_1 + k] * B_tempt[k * col_2 + j];
			}
		}
	}

	for (int i = 0; i < row_1; i++)
	{
		for (int j = 0; j < col_2; j++)
		{
			AB[i * col_2 + j] = AB_tempt[i * col_2 + j];
		}
	}

	delete[] A_tempt;
	delete[] B_tempt;
	delete[] AB_tempt;

}


// 求 [start, end] 之间的 素数， 返回值 count 表示 素数总个数
int primerange(int start, int end, int* result, int count)
{
	int* result_temp = new int[1000];
	int result_count = 0;
	bool prime;  //定义bool变量


	for (int i = start; i < (end + 1); i++)
	{
		prime = true;    //先令 prime 为真

		if ((i == 0) || (i == 1))
		{
			// cout << "i=0 or i=1" <<endl;
			continue;
		}

		for (int j = 2; j < i; j++)    //对 2 到 m 进行循环
		{
			if (i % j == 0)    //若 i 整除 j 为 0，令 prime 为假，循环终止
			{
				prime = false;
				break;
			}
		}
		if (prime)    //若 prime 为真，输出 n
		{
			result_temp[result_count] = i;
			result_count += 1;
		}

	}

	for (int i = 0; i < result_count; i++)
		result[i] = result_temp[i];

	count = result_count;

	delete[] result_temp;

	return count;

}



void gen_good_point_set(int num_particles, int dim, double* low_bound, double* up_bound, double* good_point_set)
{
	// num_particles: 群体个数； dim: 每个个体的维数； low_bound: 各个维度的 下限值；up_bound: 各个维度的 上限值； good_point_set: 生成的 佳点集

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


	for (int i = 0; i<num_particles; i++)
	{
		temp_1[i] = i + 1;

	}
	for (int i = 0; i<dim; i++)
	{
		ones_dim[i] = 1;
	}

	MatrixMultiMatrix(num_particles, 1, dim, temp_1, ones_dim, temp_1_multi_ones_dim);

	// ind = np.arange(1, dim + 1)
	for (int i = 0; i<dim; i++)
		ind[i] = i + 1;

	// prime = list(sympy.sieve.primerange(0, 100 * dim))
	prime_count = primerange(0, 100 * dim, prime_list_tempt, prime_count);

	int* prime_list = new int[prime_count];    //  注意，这是定义 prime_list，不能提前定义！
	for (int i = 0; i<prime_count; i++)
		prime_list[i] = prime_list_tempt[i];

	for (int i = 0; i<prime_count; i++)
	{
		if (prime_list[i] < (2 * dim + 3))
			continue;

		else
		{
			// idx = i;
			prime_idx = prime_list[i];
			break;
		}

	}

	for (int i = 0; i<dim; i++)
	{
		temp_2[i] = 2 * M_PI * ind[i] / prime_idx;

		temp_2[i] = 2 * cos(temp_2[i]);

	}


	for (int i = 0; i < num_particles; i++)
	{
		ones_num_particles[i] = 1;
	}

	MatrixMultiMatrix(num_particles, 1, dim, ones_num_particles, temp_2, ones_num_particles_multi_temp_2);

	//  gd = temp1 * temp2
	for (int i = 0; i < (num_particles * dim); i++)
		good_points[i] = temp_1_multi_ones_dim[i] * ones_num_particles_multi_temp_2[i];

	for (int i = 0; i < (num_particles * dim); i++)
		good_points_temp[i] = good_points[i];

	// gd = np.mod(gd, 1)
	for (int i = 0; i < (num_particles * dim); i++)
	{
		good_points_temp[i] = fmod(good_points_temp[i], 1);

		if (good_points_temp[i] < 0)
			good_points_temp[i] = 1 + good_points_temp[i];

	}

	for (int i = 0; i < (num_particles * dim); i++)
		good_points[i] = good_points_temp[i];

	/*
	aa = list(map(lambda x, y: x - y, up_bound, low_bound))
	*/
	for (int i = 0; i < dim; i++)
	{
		up_bound_minus_low_bound[i] = up_bound[i] - low_bound[i];
	}


	for (int i = 0; i < num_particles; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			//  bb = list(map(lambda x, y: x * y, [gd[:, i] for i in range(dim)], aa))
			good_points[i*dim + j] = good_points[i*dim + j] * up_bound_minus_low_bound[j];

			//  cc = list(map(lambda x, y: x + y, bb, low_bound))
			good_points[i*dim + j] = good_points[i*dim + j] + low_bound[j];
		}
	}

	for (int i = 0; i < (num_particles * dim); i++)
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


