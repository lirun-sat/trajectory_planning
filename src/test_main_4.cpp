#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include"data.h"
#include "utils.h"
#include "WOAAlgorithm.h"
#include "TLBO_Algorithm.h"


double calc_fitness_woa(WOA_Whale& whale);
double calc_fitness_tlbo(Student& student);
void sort_max2min_main(double*, int, int*);


// 仅仅考虑 末端 最终位置 和 基座姿态变化， 没有考虑 执行时间 和 路径是否是直线
int main()
{
	clock_t start;
	clock_t end;

	cout << "Initial processing..." << endl;
	cout << "delta_tau = " << "  " << delta_tau << endl;
	cout << "RPY_END_DESIRED:" << endl;
	cout << RPY_END_DESIRED[0] << "  " << RPY_END_DESIRED[1] << "  " << RPY_END_DESIRED[2] << endl;
	cout << "Pe_DESIRED:" << endl;
	cout << Pe_DESIRED[0] << "  " << Pe_DESIRED[1] << "  " << Pe_DESIRED[2] << endl;


	/*注意，徐 论文中的 3-41 和 3-42 有问题
	末端初始位置 根据 初始关节角和原点定义（整个系统质心）可以求出，末端的初始姿态 为 【0.0000， 0.9848， 0.0000， 0.1736】
	末端的期望位置 [1.3501, -0.0137, 0.1401], 末端的期望姿态 [-0.0085, 0.9988, 0.0096, 0.0472]
    */

	int dimension = N;
	int student_whale_count = 40;
	// int Iter_Max = 40000;
	int Iter_Max = 100000;
	int num_calc = 1000;
	
	cout << "Iter_Max = " << "  " << Iter_Max << endl;
	cout << "num_calc = " << "  " << num_calc << endl;
	cout << "student_whale_count = " << "  " << student_whale_count << endl;

	double* result_min = new double[N];
	double* result_max = new double[N];
	double* para_low_bound = new double[N];
	double* para_up_bound = new double[N];
	double* woa_tlbo_fitness = new double[student_whale_count];
	int* woa_tlbo_index = new int[student_whale_count];
	
	double** woa_tlbo_position = new double*[student_whale_count];
	for (int i = 0; i < student_whale_count; i++)
	{
		woa_tlbo_position[i] = new double[dimension];
	}
	
	double woa_tlbo_fitness_temp_1 = 0;
	double woa_tlbo_fitness_temp_2 = 0;
	double woa_tlbo_fitness_temp_3 = 0;
	double woa_tlbo_fitness_temp_4 = 0;
	double woa_tlbo_fitness_temp_5 = 0;
	double woa_tlbo_fitness_temp_6 = 0;
	double woa_tlbo_fitness_temp_7 = 0;
	double woa_tlbo_fitness_temp_8 = 0;
	double woa_tlbo_fitness_temp_9 = 0;
	double woa_tlbo_fitness_temp_15 = 0;
	double woa_tlbo_fitness_temp_22 = 0;
	double woa_tlbo_fitness_temp_35 = 0;
	double woa_tlbo_fitness_temp_50 = 0;
	double woa_tlbo_fitness_temp_80 = 0;

	
// *******************************************************************************************************************************************************************
	calc_para_range(result_min, result_max);

	for (int i = 0; i<N; i++)
	{
		para_low_bound[i] = result_min[i];
		para_up_bound[i] = result_max[i];
	}

	cout << "result_min:" << endl;
	for(int i=0;i<N;i++)
	{
		cout << result_min[i] << "  ";
	}
	cout << endl;

	cout << "result_max:" << endl;
	for(int i=0;i<N;i++)
	{
		cout << result_max[i] << "  ";
	}
	cout << endl;
// *******************************************************************************************************************************************************************

	//构建算法
	WOA_Algorithm woa(calc_fitness_woa, para_low_bound, para_up_bound, dimension, student_whale_count / 2);
	woa.randomlyInitial();
	woa.refresh();

	TLBO_Algorithm tlbo(calc_fitness_tlbo, para_low_bound, para_up_bound, dimension, student_whale_count / 2);
	tlbo.randomlyInitial();
	tlbo.refresh();

// *****************************************************************************************************************************************************************
	cout << "算法构建完成，进入迭代搜索过程......" << endl;

	for (int kk = 0; kk < num_calc; kk++)
	{	
		cout << "第" << kk << "次搜索......" << endl;
		start = clock();

		woa.randomlyInitial();
		cout << "WOA 初始化完成......" << endl;
		woa.refresh();
		tlbo.randomlyInitial();
		cout << "TLBO 初始化完成......" << endl;
		tlbo.refresh();

		for (int Iter = 0; Iter < Iter_Max; Iter++)
		{
			for (int i = 0; i < student_whale_count; i++)
			{
				if (i < student_whale_count / 2)
				{
					woa_tlbo_fitness[i] = woa._whaleSet[i]._fitness;
					for (int j = 0; j < dimension; j++)
					{
						woa_tlbo_position[i][j] = woa._whaleSet[i]._position[j];
					}
				}
				else
				{
					woa_tlbo_fitness[i] = tlbo._studentSet[i - student_whale_count / 2]._fitness;
					for (int j = 0; j < dimension; j++)
					{
						woa_tlbo_position[i][j] = tlbo._studentSet[i - student_whale_count / 2]._position[j];
					}
				}
				woa_tlbo_index[i] = i;
			}


// ***********************************************  Termination condition  *************************************************************************    
			sort_max2min_main(woa_tlbo_fitness, student_whale_count, woa_tlbo_index);
			
			if (Iter == Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;	
				woa_tlbo_fitness_temp_1 = woa_tlbo_fitness[student_whale_count - 1];
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_1 << endl;
				if(fabs(woa_tlbo_fitness_temp_1 - woa_tlbo_fitness_temp_2) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
			}
			else if (Iter == 2*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;	
				woa_tlbo_fitness_temp_2 = woa_tlbo_fitness[student_whale_count - 1];
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_2 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_2 - woa_tlbo_fitness_temp_1) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			else if (Iter == 3*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;	
				woa_tlbo_fitness_temp_3 = woa_tlbo_fitness[student_whale_count - 1];
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_3 << endl;
				if(woa_tlbo_fitness_temp_3 > 200)
				{
					cout << "Too large woa_tlbo_fitness after 300 steps " << endl;
					break;
				}
				/*
				if(fabs(woa_tlbo_fitness_temp_3 - woa_tlbo_fitness_temp_2) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				*/
			}
			/*
			else if (Iter == 4*Iter_Max / 40)
			{
				cout << "Iteration is :" << Iter << endl;	
				woa_tlbo_fitness_temp_4 = woa_tlbo_fitness[student_whale_count - 1];
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_4 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_4 - woa_tlbo_fitness_temp_3) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			else if (Iter == 5*Iter_Max / 40)
			{
				cout << "Iteration is :" << Iter << endl;	
				woa_tlbo_fitness_temp_5 = woa_tlbo_fitness[student_whale_count - 1];
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_5 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_5 - woa_tlbo_fitness_temp_4) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			*/
			else if (Iter == 6*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;	
				woa_tlbo_fitness_temp_6 = woa_tlbo_fitness[student_whale_count - 1];
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_6 << endl;
				
				// if(fabs(woa_tlbo_fitness_temp_6 - woa_tlbo_fitness_temp_5) < 0.01)
				if(fabs(woa_tlbo_fitness_temp_6 - woa_tlbo_fitness_temp_3) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			/*
			else if (Iter == 7*Iter_Max / 40)
			{
				cout << "Iteration is :" << Iter << endl;	
				woa_tlbo_fitness_temp_7 = woa_tlbo_fitness[student_whale_count - 1];
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_7 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_7 - woa_tlbo_fitness_temp_6) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			else if (Iter == 8*Iter_Max / 80)
			{
				cout << "Iteration is :" << Iter << endl;	
				woa_tlbo_fitness_temp_8 = woa_tlbo_fitness[student_whale_count - 1];
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_8 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_8 - woa_tlbo_fitness_temp_7) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			*/
			else if (Iter == 9*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;
				woa_tlbo_fitness_temp_9 = woa_tlbo_fitness[student_whale_count - 1];	
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_9 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_9 - woa_tlbo_fitness_temp_6) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			
			else if (Iter == 15*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;
				woa_tlbo_fitness_temp_15 = woa_tlbo_fitness[student_whale_count - 1];	
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_15 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_15 - woa_tlbo_fitness_temp_9) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			
			else if (Iter == 22*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;
				woa_tlbo_fitness_temp_22 = woa_tlbo_fitness[student_whale_count - 1];	
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_22 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_22 - woa_tlbo_fitness_temp_15) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			
			else if (Iter == 35*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;
				woa_tlbo_fitness_temp_35 = woa_tlbo_fitness[student_whale_count - 1];	
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_35 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_35 - woa_tlbo_fitness_temp_22) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			else if (Iter == 50*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;
				woa_tlbo_fitness_temp_50 = woa_tlbo_fitness[student_whale_count - 1];	
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_50 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_50 - woa_tlbo_fitness_temp_35) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			else if (Iter == 80*Iter_Max / 100)
			{
				cout << "Iteration is :" << Iter << endl;
				woa_tlbo_fitness_temp_80 = woa_tlbo_fitness[student_whale_count - 1];	
				cout << "woa_tlbo_fitness is :" << woa_tlbo_fitness_temp_80 << endl;
				
				if(fabs(woa_tlbo_fitness_temp_80 - woa_tlbo_fitness_temp_50) < 0.01)
				{
					cout << "Convergence reached" << endl;
					break;
				}
				
			}
			if (woa_tlbo_fitness[student_whale_count - 1] < 1)
			{
				cout << "Find solution, Iter:" << Iter << endl;
				break;
			}

//***************************************************************************************************************************************************

			for (int i = 0; i < student_whale_count; i++)
			{
				if (i < student_whale_count / 2)
				{
					woa._whaleSet[i]._fitness = woa_tlbo_fitness[i];

					for (int j = 0; j < dimension; j++)
					{
						woa._whaleSet[i]._position[j] = woa_tlbo_position[woa_tlbo_index[i]][j];
					}

				}
				else
				{
					tlbo._studentSet[i - student_whale_count / 2]._fitness = woa_tlbo_fitness[i];

					for (int j = 0; j < dimension; j++)
					{
						tlbo._studentSet[i - student_whale_count / 2]._position[j] = woa_tlbo_position[woa_tlbo_index[i]][j];
					}
				}

			}

			woa.refresh();
			tlbo.refresh();
			woa.update();
			tlbo.update();

		}

		for (int i = 0; i < student_whale_count; i++)
		{
			if (i < student_whale_count / 2)
			{
				woa_tlbo_fitness[i] = woa._whaleSet[i]._fitness;
				for (int j = 0; j < dimension; j++)
				{
					woa_tlbo_position[i][j] = woa._whaleSet[i]._position[j];
				}
			}
			else
			{
				woa_tlbo_fitness[i] = tlbo._studentSet[i - student_whale_count / 2]._fitness;
				for (int j = 0; j < dimension; j++)
				{
					woa_tlbo_position[i][j] = tlbo._studentSet[i - student_whale_count / 2]._position[j];
				}
			}
			
			woa_tlbo_index[i] = i;

		}

		// for (int k = 0; k < student_whale_count; k++)
		// 	woa_tlbo_index[k] = k;

		sort_max2min_main(woa_tlbo_fitness, student_whale_count, woa_tlbo_index);

		cout << "woa_tlbo_fitness:" << "  " << woa_tlbo_fitness[student_whale_count - 1] << endl;
		for (int j = 0; j < dimension; j++)
		{
			cout << woa_tlbo_position[woa_tlbo_index[student_whale_count - 1]][j] << "  ";
		}
		cout << endl;

		end = clock();   //结束时间
		cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;  //输出时间（单位：ｓ）

	}


	delete[] woa_tlbo_fitness;
	delete[] woa_tlbo_index;

	for (int i = 0; i<student_whale_count; i++)
		delete[] woa_tlbo_position[i];
	delete[] woa_tlbo_position;

	delete[] para_low_bound;
	delete[] para_up_bound;
	delete[] result_min;
	delete[] result_max;

}



// cost function 只考虑末端位姿
double calc_fitness_woa(WOA_Whale& whale)
{
	double* para = new double[N];
	double* eta_end = new double;
	double* xi_end = new double[3];
	double* Pe = new double[3];
	double* eta_b = new double;
	double* xi_b = new double[3];
	double* quaternion_end_desired = new double[4];
	double* eta_end_desired = new double;
	double* xi_end_desired = new double[3];
	double* delta_eta_base = new double;
	double* delta_eta_end = new double;
	double* delta_xi_base = new double[3];
	double* delta_xi_end = new double[3];
	double* delta_Pe_end = new double[3];
    double* p_e_initial = new double[3];
    double* locus = new double;
    double* delta_xi_b_distrb_max = new double;
    double* manipl = new double;
    double* T_min = new double;

	for (int i = 0; i < N; i++)
		para[i] = whale._position[i];

	// double K_a = 1 / 0.0087;
    double K_a = 1 / 0.0002;  // 1 degree error tolerance
	// double K_p = 1 / 0.005;  //  0.005 meter error tolerance
	double K_p = 1 / 0.002;  //  0.002 meter error tolerance for end-effector
    // double K_s = 1 / 0.1;  // locus tolerance 0.1 meter.
    double K_s = 0;  // locus tolerance 0.1 meter.
    double K_b = 1 / 0.0008;  // attitude error tolerance for base is 5 degree
    double K_M = 0;
    double K_t = 0;  //  max allowed motion time is set to 100 seconds, tolerance is 10 seconds.

	double cost_func = 0;

	double delta_xi_end_mod_temp = 0;
	double delta_Pe_end_mod_temp = 0;
	double delta_xi_base_mod_temp = 0;

	forward_kin_2( para,  eta_end,  xi_end,  Pe,  eta_b,  xi_b,  p_e_initial,  locus,  delta_xi_b_distrb_max,  manipl, T_min);

// ********************************************************  straight_line_locus  ************************************************************************************
    double delta_p_e = 0;
    double straight_line_locus = 0;
    for(int i=0;i<3;i++)
    {
        delta_p_e = p_e_initial[i] - Pe[i];
        straight_line_locus += delta_p_e * delta_p_e;
    }
    straight_line_locus = sqrt(straight_line_locus);
// **********************************************************************************************************************************************************************

	zyx2quaternion(RPY_END_DESIRED[2], RPY_END_DESIRED[1], RPY_END_DESIRED[0], quaternion_end_desired);

	*eta_end_desired = quaternion_end_desired[0];

	for (int i = 0; i<3; i++)
	{
		xi_end_desired[i] = quaternion_end_desired[i + 1];
	}

	delta_var(Pe_DESIRED, eta_end_desired, xi_end_desired, Pe, *eta_end, xi_end, *eta_b, xi_b, delta_eta_end, delta_xi_end, delta_Pe_end, delta_eta_base, delta_xi_base);

	for (int i = 0; i < 3; i++)
	{
		delta_xi_end_mod_temp += delta_xi_end[i] * delta_xi_end[i];
		delta_Pe_end_mod_temp += delta_Pe_end[i] * delta_Pe_end[i];
		delta_xi_base_mod_temp += delta_xi_base[i] * delta_xi_base[i];
	}

	delta_xi_end_mod_temp = sqrt(delta_xi_end_mod_temp); 
    delta_Pe_end_mod_temp = sqrt(delta_Pe_end_mod_temp); 
    delta_xi_base_mod_temp = sqrt(delta_xi_base_mod_temp);


	// ****************************************  RPY error of end-effector,  Pe error  ********************************************************************************

	// cost_func = K_a * delta_xi_end_mod_temp + K_p * delta_Pe_end_mod_temp;

    cost_func = K_a * delta_xi_end_mod_temp 
                + K_p * delta_Pe_end_mod_temp 
                + K_b * (*delta_xi_b_distrb_max) 
                + K_s * fabs((*locus) - straight_line_locus)  
                + K_M * (*manipl) 
                + K_t * (*T_min);


    delete[] para;
    delete eta_end;
    delete[] xi_end ;
    delete[] Pe ;
    delete eta_b ;
	delete[] xi_b ;
	delete[] quaternion_end_desired;
	delete eta_end_desired;
	delete[] xi_end_desired ;
	delete delta_eta_base ;
	delete delta_eta_end ;
	delete[] delta_xi_base ;
	delete[] delta_xi_end ;
	delete[] delta_Pe_end ; 
	delete[] p_e_initial ;
	delete locus ;
	delete delta_xi_b_distrb_max;
	delete manipl;
	delete T_min;


	return cost_func;

}




double calc_fitness_tlbo(Student& student)
{
	double* para = new double[N];
	double* eta_end = new double;
	double* xi_end = new double[3];
	double* Pe = new double[3];
	double* eta_b = new double;
	double* xi_b = new double[3];
	double* quaternion_end_desired = new double[4];
	double* eta_end_desired = new double;
	double* xi_end_desired = new double[3];
	double* delta_eta_base = new double;
	double* delta_eta_end = new double;
	double* delta_xi_base = new double[3];
	double* delta_xi_end = new double[3];
	double* delta_Pe_end = new double[3];
    double* p_e_initial = new double[3];
    double* locus = new double;
    double* delta_xi_b_distrb_max = new double;
    double* manipl = new double;
    double* T_min = new double;

	for (int i = 0; i < N; i++)
		para[i] = student._position[i];

	// double K_a = 1 / 0.0087;
    double K_a = 1 / 0.0002;  // 1 degree error tolerance for end-effector
	// double K_p = 1 / 0.005;  //  0.005 meter error tolerance for end-effector
	double K_p = 1 / 0.002;  //  0.002 meter error tolerance for end-effector
    // double K_s = 1 / 0.1;  // locus tolerance 0.1 meter for end-effector.
    double K_s = 0;  // locus tolerance 0.1 meter for end-effector.
    double K_b = 1 / 0.0008;  // attitude error tolerance for base is 5 degree
    double K_M = 0;
    double K_t = 0;  //  max allowed motion time is set to 100 seconds, tolerance is 10 seconds.

	double cost_func = 0;

	double delta_xi_end_mod_temp = 0;
	double delta_Pe_end_mod_temp = 0;
	double delta_xi_base_mod_temp = 0;

	forward_kin_2( para,  eta_end,  xi_end,  Pe,  eta_b,  xi_b,  p_e_initial,  locus,  delta_xi_b_distrb_max,  manipl, T_min);

// ********************************************************  straight_line_locus  ************************************************************************************
    double delta_p_e = 0;
    double straight_line_locus = 0;
    for(int i=0;i<3;i++)
    {
        delta_p_e = p_e_initial[i] - Pe[i];
        straight_line_locus += delta_p_e * delta_p_e;
    }
    straight_line_locus = sqrt(straight_line_locus);
// **********************************************************************************************************************************************************************

	zyx2quaternion(RPY_END_DESIRED[2], RPY_END_DESIRED[1], RPY_END_DESIRED[0], quaternion_end_desired);

	*eta_end_desired = quaternion_end_desired[0];

	for (int i = 0; i<3; i++)
	{
		xi_end_desired[i] = quaternion_end_desired[i + 1];
	}

	delta_var(Pe_DESIRED, eta_end_desired, xi_end_desired, Pe, *eta_end, xi_end, *eta_b, xi_b, delta_eta_end, delta_xi_end, delta_Pe_end, delta_eta_base, delta_xi_base);

	for (int i = 0; i < 3; i++)
	{
		delta_xi_end_mod_temp += delta_xi_end[i] * delta_xi_end[i];
		delta_Pe_end_mod_temp += delta_Pe_end[i] * delta_Pe_end[i];
		delta_xi_base_mod_temp += delta_xi_base[i] * delta_xi_base[i];
	}

	delta_xi_end_mod_temp = sqrt(delta_xi_end_mod_temp);
	delta_Pe_end_mod_temp = sqrt(delta_Pe_end_mod_temp);
	delta_xi_base_mod_temp = sqrt(delta_xi_base_mod_temp);


    cost_func = K_a * delta_xi_end_mod_temp 
                + K_p * delta_Pe_end_mod_temp 
                + K_b * (*delta_xi_b_distrb_max) 
                + K_s * fabs(*locus - straight_line_locus)  
                + K_M * (*manipl) 
                + K_t * (*T_min);



    delete[] para ;
    delete eta_end ;
    delete[] xi_end ;
    delete[] Pe ;
    delete eta_b ;
	delete[] xi_b ;
	delete[] quaternion_end_desired; 
	delete eta_end_desired ;
	delete[] xi_end_desired ;
	delete delta_eta_base ;
	delete delta_eta_end ;
	delete[] delta_xi_base ;
	delete[] delta_xi_end ;
	delete[] delta_Pe_end ;
	delete[] p_e_initial;
	delete locus ;
	delete delta_xi_b_distrb_max ;
	delete manipl;
	delete T_min ;


	return cost_func;

}



// 数组排序，从大到小排序，并返回排序后的数组对应原数组的下标
void sort_max2min_main(double* a, int length, int* b)
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














