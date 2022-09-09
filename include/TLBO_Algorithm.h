
#ifndef _TLBOALGORITHM_H
#define _TLBOALGORITHM_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iomanip>

#define EPS 0.000001


// double avg(double* parameter, int n);
// double stddev(double* parameter, int n);
void MatrixMultiMatrix(int row_1, int col_1, int col_2, double* A, double* B, double* AB);
int primerange(int start, int end, int* result, int count);
void gen_good_point_set(int num_particles, int dim, double* low_bound, double* up_bound, double* good_point_set);
void bubbleSort(double* arr, int n);



//粒子群算法例子个体
class Student
{
public:
	int _dimension;
	double* _position;
	double _fitness;

	//构造函数，粒子维度初始化为0
	Student(void){
		_dimension = 0;
	}

	//析构函数，释放粒子内存
	~Student(void){
		if (_dimension){
			delete[] _position;
		}
	}

	//初始化函数，用于为粒子开辟内存空间
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
	//复制函数，用于粒子间的复制操作
	void copy(Student& student){
		this->initial(student._dimension);
		for (int i = 0; i < _dimension; ++i){
			_position[i] = student._position[i];
		}
		_fitness = student._fitness;

	}

};


//PSO算法
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
	double (*_fitnessFunction)(Student&);

	TLBO_Algorithm(double(*objFunction)(Student&), double *positionMinValue, double *positionMaxValue, int dimension, int studentCount){
		//初始化类内参数并分配内存
		_fitnessFunction = objFunction;
		_dimension = dimension;
		_studentCount = studentCount;
		_positionMinValue = new double[_dimension];
		_positionMaxValue = new double[_dimension];
		_studentSet = new Student[_studentCount];

		for (int i = 0; i < _dimension; ++i){
			_positionMinValue[i] = positionMinValue[i];
			_positionMaxValue[i] = positionMaxValue[i];
		}
		for (int i = 0; i < _studentCount; ++i){
			_studentSet[i].initial(_dimension);
		}
		_best_student.initial(_dimension);

		//配置随机数种子
		srand((unsigned int)time(NULL));

	}
	/***************************************************************
	* 函数描述：析构，释放算法内存
	***************************************************************/
	~TLBO_Algorithm(void){
		//释放内存
		delete[] _positionMinValue;
		delete[] _positionMaxValue;
		delete[] _studentSet;

	}

	double rand0_1(void){
		return((1.0 * rand()) / RAND_MAX);
	}

	void randomlyInitial(void){
		// double* good_point_set_temp;
		// good_point_set_temp = new double[_studentCount * _dimension];
		//  产生佳点集
		// gen_good_point_set(_studentCount, _dimension, _positionMinValue, _positionMaxValue, good_point_set_temp);

		for (int i = 0; i < _studentCount; ++i){
			_studentSet[i]._position[0] = 0.852354 + 0.1 * rand0_1();
			_studentSet[i]._position[1] = 0.433124 + 0.1 * rand0_1();
			_studentSet[i]._position[2] = 2.03323  + 0.5 * rand0_1();
			_studentSet[i]._position[3] = 2.9251   + 0.5 * rand0_1();
			_studentSet[i]._position[4] = 1.82466  + 0.5 * rand0_1();
			_studentSet[i]._position[5] = 1.94137  + 0.5 * rand0_1();
			_studentSet[i]._position[6] = 0.966595 + 0.5 * rand0_1();
			
			_studentSet[i]._fitness = _fitnessFunction(_studentSet[i]);
		}

		// delete[] good_point_set_temp;

	}

	// 选出最小的适应度函数值，找到对应的粒子，将其作为最优粒子
	void refresh(void){
		_globalBestStudentIndex = -1;
		double* fitness_temp = new double[_studentCount];

		for (int i = 0; i < _studentCount; ++i){
			fitness_temp[i] = _studentSet[i]._fitness;
		}

		bubbleSort(fitness_temp, _studentCount);

		// 选出适应度函数最小的值
		_globalBestFitness = fitness_temp[_studentCount - 1];

		for (int i = 0; i < _studentCount; ++i){
			if ((_studentSet[i]._fitness - _globalBestFitness) < EPS){
				_globalBestStudentIndex = i;
			}
			else{
				continue;
			}
		}
		if (_globalBestStudentIndex != -1){
			_best_student.copy(_studentSet[_globalBestStudentIndex]);
		}

		delete[] fitness_temp;

	}


	void update(void){
		double* X_mean = new double[_dimension];
		double X_mean_temp;

		Student Xnew_temp_1;
		Xnew_temp_1.initial(_dimension);

		Student Xnew_temp_2;
		Xnew_temp_2.initial(_dimension);

		for (int i = 0; i < _studentCount; ++i){
			for (int j = 0; j < _dimension; ++j){
				X_mean_temp = 0;
				for (int k = 0; k < _studentCount; k++){
					X_mean_temp += _studentSet[k]._position[j];
				}
				X_mean[j] = X_mean_temp / _studentCount;
			}

			//  要取得[a,b]的随机整数，使用(rand() % (b-a+1))+ a;
			TF = (rand() % (2 - 1 + 1)) + 1;   

			// Teaching phase-----------
			for (int j = 0; j < _dimension; ++j){
				r1 = rand0_1();
				Xnew_temp_1._position[j] = _studentSet[i]._position[j] + r1 * (_best_student._position[j] - TF * X_mean[j]);
				if (Xnew_temp_1._position[j] >= _positionMaxValue[j]){
					Xnew_temp_1._position[j] = (_studentSet[i]._position[j] + _positionMaxValue[j]) / 2;
				}
				else if (Xnew_temp_1._position[j] <= _positionMinValue[j]){
					Xnew_temp_1._position[j] = (_studentSet[i]._position[j] + _positionMinValue[j]) / 2;
				}
			}
			Xnew_temp_1._fitness = this->_fitnessFunction(Xnew_temp_1);
			if (Xnew_temp_1._fitness < _studentSet[i]._fitness){
				_studentSet[i].copy(Xnew_temp_1);
			}
			if (Xnew_temp_1._fitness < _best_student._fitness){
				_best_student.copy(Xnew_temp_1);
			}

			// Learning phase----------------------------------------------------------------------------------------------------------------------------------------------
			q = (rand() % ((_studentCount - 1) - 0 + 1)) + 0;   //  要取得[a,b]的随机整数，使用(rand() % (b-a+1))+ a;

			while (q == i){
				q = (rand() % ((_studentCount - 1) - 0 + 1)) + 0;
			}
			if (_studentSet[i]._fitness < _studentSet[q]._fitness){
				for (int j = 0; j < _dimension; ++j){
					r2 = rand0_1();
					Xnew_temp_2._position[j] = _studentSet[i]._position[j] + r2 * (_studentSet[i]._position[j] - _studentSet[q]._position[j]);
				}
			}
			else{
				for (int j = 0; j < _dimension; ++j){
					r2 = rand0_1();
					Xnew_temp_2._position[j] = _studentSet[i]._position[j] - r2 * (_studentSet[i]._position[j] - _studentSet[q]._position[j]);
				}
			}
			for (int j = 0; j < _dimension; ++j){
				if (Xnew_temp_2._position[j] >= _positionMaxValue[j]){
					Xnew_temp_2._position[j] = (_studentSet[i]._position[j] + _positionMaxValue[j] / 2);
				}
				else if (Xnew_temp_2._position[j] <= _positionMinValue[j]){
					Xnew_temp_2._position[j] = (_studentSet[i]._position[j] + _positionMinValue[j] / 2);
				}
			}

			Xnew_temp_2._fitness = this->_fitnessFunction(Xnew_temp_2);

			if (Xnew_temp_2._fitness < _studentSet[i]._fitness){
				_studentSet[i].copy(Xnew_temp_2);
			}
			if (Xnew_temp_2._fitness < _best_student._fitness){
				_best_student.copy(Xnew_temp_2);
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





#endif 


