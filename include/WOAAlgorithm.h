#ifndef _WOAALGORITHM_H
#define _WOAALGORITHM_H

#define _USE_MATH_DEFINES 

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iomanip>
#include <random>
#include <chrono>


void MatrixMultiMatrix(int row_1, int col_1, int col_2, double* A, double* B, double* AB);
int primerange(int start, int end, int* result, int count);
void gen_good_point_set(int num_particles, int dim, double* low_bound, double* up_bound, double* good_point_set);
void bubbleSort(double* arr, int n);


class WOA_Whale
{
public:
	int _dimension;     //维度
	double* _position;  //所在位置数组指针
	double _fitness;    //适应度

	//构造函数，维度初始化为0
	WOA_Whale(void){
		_dimension = 0;
	}
	//析构函数，释放内存
	~WOA_Whale(void){
		if (_dimension){
			delete[] _position;
		}
	}
	//初始化函数，用于开辟内存空间
	void initial(int dimension){
		if (_dimension != dimension && dimension){
			//需要重新分配内存
			if (_dimension){
				//消除已有内存
				delete[] _position;
			}
			//开辟新内存
			_dimension = dimension;
			_position = new double[_dimension];
		}
	}
	//复制函数，用于复制操作
	void copy(WOA_Whale& whale){
		this->initial(whale._dimension);
		for (int i = 0; i<_dimension; ++i){
			_position[i] = whale._position[i];
			_fitness = whale._fitness;
		}
	}
};

//算法
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
	double (*_fitnessFunction)(WOA_Whale&);

	WOA_Algorithm(double(*objFunction)(WOA_Whale&), double *positionMinValue, double *positionMaxValue, int dimension, int whaleCount){
		//初始化类内参数并分配内存
		_fitnessFunction = objFunction;
		_dimension = dimension;
		_positionMinValue = new double[_dimension];
		_positionMaxValue = new double[_dimension];
		_whaleCount = whaleCount;
		_whaleSet = new WOA_Whale[_whaleCount];

		for (int i = 0; i < _dimension; ++i){
			_positionMinValue[i] = positionMinValue[i];
			_positionMaxValue[i] = positionMaxValue[i];
		}
		for (int i = 0; i < _whaleCount; ++i){
			_whaleSet[i].initial(_dimension);
		}
		_best_Whale.initial(_dimension);
		//配置随机数种子
		srand((unsigned int)time(NULL));
	}
	/***************************************************************
	* 函数描述：析构一个PSO算法，释放算法内存
	***************************************************************/
	~WOA_Algorithm(void){
		//释放内存
		delete[] _positionMinValue;
		delete[] _positionMaxValue;
		delete[] _whaleSet;
	}

	double rand0_1(void){
		return((1.0 * rand()) / RAND_MAX);
	}

	double normal_distribution_number(double mean_temp, double std_devation){
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator(seed);
		std::normal_distribution<double> distribution(mean_temp, std_devation);

		return(distribution(generator));
	}

	double levy_flight_step(){
		// double beta = rand0_1() * 2;  // beta is stochastic 
		double beta = 1.5;
		double alpha_u = pow(((tgamma(1+beta)*sin(M_PI*beta/2)) / (tgamma((1+beta)/2) * beta * pow(2, (beta-1)/2))), (1/beta));
		double alpha_v = 1;
		double u = normal_distribution_number(0, alpha_u);
		double v = normal_distribution_number(0, alpha_v);
		double step = u / pow(fabs(v), (1 / beta));

		return step;

	}

	double sign_fun(double x){
		if(x==0){
			return 0;
		}
		else if(x>0){
			return 1;
		}
		else{
			return (-1);
		}
	}

	void randomlyInitial(void){
		// double* good_point_set_temp;
		// good_point_set_temp = new double[_whaleCount* _dimension];
		// //  产生佳点集
		// gen_good_point_set(_whaleCount, _dimension, _positionMinValue, _positionMaxValue, good_point_set_temp);
		// for (int i = 0; i < _whaleCount; ++i){
		// 	for (int j = 0; j < _dimension; ++j){
		// 		_whaleSet[i]._position[j] = good_point_set_temp[i * _dimension + j];
		// 	}
		// 	_whaleSet[i]._fitness = _fitnessFunction(_whaleSet[i]);
		// }
		// delete[] good_point_set_temp;

		for (int i = 0; i < _whaleCount; ++i){
			for (int j = 0; j < _dimension; ++j){
				_whaleSet[i]._position[j] = _positionMinValue[j] + rand0_1() * (_positionMaxValue[j] - _positionMinValue[j]);
			}
			_whaleSet[i]._fitness = _fitnessFunction(_whaleSet[i]);
		}

	}

	void refresh(void){
		_globalBestWhaleIndex = -1;
		double* fitness_temp = new double[_whaleCount];

		for (int i = 0; i < _whaleCount; ++i){
			fitness_temp[i] = _whaleSet[i]._fitness;
		}
		// fitness_temp 从大到小排序, sorting is not necessary
		bubbleSort(fitness_temp, _whaleCount);  
		// 选出适应度函数最小的值
		_globalBestFitness = fitness_temp[_whaleCount - 1];
		for (int i = 0; i < _whaleCount; ++i){
			if (_whaleSet[i]._fitness <= _globalBestFitness){
				_globalBestWhaleIndex = i;
			}
			else{
				continue;
			}
		}
		if (_globalBestWhaleIndex != -1){
			_best_Whale.copy(_whaleSet[_globalBestWhaleIndex]);
		}
			
		delete[] fitness_temp;
	}

	void update(int iter, int max_iter){
		a = 2 * (1 - (double)iter / (double)max_iter);  
		//  需要对 iter 进行类型转换， 因为 a 是 double 型变量.
		a2 = (-1) + (double)iter * ((-1) / (double)max_iter);

		double* X_D = new double[_dimension];
		double* Xnew = new double[_dimension];
		double* X_D_prime = new double[_dimension];
		double* Xrand = new double[_dimension];
		double r1;
		double r2;
		double r3;
		double p;
		int q;

		for (int i = 0; i < _whaleCount; ++i){
			r1 = rand0_1();
			r2 = rand0_1();
			r3 = rand0_1();
			b = 1;
			A = a * (2 * r1 - 1);
			C = 2 * r2;
			l = (a2 - 1) * r3 + 1;
			p = rand0_1();

			if (p < 0.5){
				if (fabs(A) > 1){
					q = (rand() % ((_whaleCount - 1) - 0 + 1)) + 0;   
					//  要取得[a,b]的随机整数，使用(rand() % (b-a+1))+ a;
					// q = (rand() % (_whaleCount - 0)) + 0;  // 要取得[a,b)的随机整数，使用(rand() % (b-a))+ a;
					while (q == i){
						q = (rand() % ((_whaleCount - 1) - 0 + 1)) + 0;
					}	
					for (int j = 0; j < _dimension; ++j){
						Xrand[j] = _whaleSet[q]._position[j];
						X_D[j] = fabs(C * Xrand[j] - _whaleSet[i]._position[j]);
						Xnew[j] = Xrand[j] - A * X_D[j];
						Xnew[j] = Xnew[j] + rand0_1() * sign_fun(rand0_1() - 0.5) * levy_flight_step();   // ************************************************
					}
				}
				else{
					for (int j = 0; j < _dimension; ++j){
						X_D[j] = fabs(C * _best_Whale._position[j] - _whaleSet[i]._position[j]);
						Xnew[j] = _best_Whale._position[j] - A * X_D[j];
						Xnew[j] = Xnew[j] + rand0_1() * sign_fun(rand0_1() - 0.5) * levy_flight_step();  // **********************************************
					}
				}
			}
			else{
				for (int j = 0; j < _dimension; ++j){
					X_D_prime[j] = fabs(_best_Whale._position[j] - _whaleSet[i]._position[j]);
					Xnew[j] = X_D_prime[j] * exp(b * l) * cos(2 * M_PI * l) + _best_Whale._position[j];
					Xnew[j] = Xnew[j] + rand0_1() * sign_fun(rand0_1() - 0.5) * levy_flight_step();     // **************************************************
				}
			}
			for (int j = 0; j < _dimension; ++j){
				if (Xnew[j] <= _positionMaxValue[j] && Xnew[j] >= _positionMinValue[j]){
					_whaleSet[i]._position[j] = Xnew[j];
				}
				else if (Xnew[j] > _positionMaxValue[j]){
					_whaleSet[i]._position[j] = (_whaleSet[i]._position[j] + _positionMaxValue[j]) / 2;
				}
				else{
					_whaleSet[i]._position[j] = (_whaleSet[i]._position[j] + _positionMinValue[j]) / 2;
				}
			}
			_whaleSet[i]._fitness = this->_fitnessFunction(_whaleSet[i]);
		}

		delete[] X_D;
		delete[] Xnew;
		delete[] X_D_prime;
		delete[] Xrand;

	}

};


#endif // _ZPSOALGORITHM_H
