#define _USE_MATH_DEFINES 
#include "forward_kinematics.h"
#include "data.h"
#include "utils.h"
#include "PSOAlgorithm.h"


double calc_fitness_pso(PSO_Particle& particle);

int main()
{
	clock_t start;
	clock_t end;

	int dimension = N;
	int particleCount = 50;
	int Iter_Max = 2000;
	int num_calc = 100;

	double inertGuideCoe = 0.9;
	double globalGuideCoe = 1.5;
	double localGuideCoe = 1.5;
	double disturbanceRate = 0.2;
	double disturbanceVelocityCoe = 0.05;
	PSO_Particle bestParticle;

	double* maxSpeed = new double[N];
	double* minSpeed = new double[N];
	maxSpeed[0] = 1.95;
	maxSpeed[1] = 0.628;
	maxSpeed[2] = 1.95;
	maxSpeed[3] = 1.95;
	maxSpeed[4] = 1.95;
	maxSpeed[5] = 0.8;
	maxSpeed[6] = 1.95;

	minSpeed[0] = -1.95;
	minSpeed[1] = -0.628;
	minSpeed[2] = -1.95;
	minSpeed[3] = -1.95;
	minSpeed[4] = -1.95;
	minSpeed[5] = -0.8;
	minSpeed[6] = -1.95;

	cout << "Initial processing..." << endl;
	cout << "delta_tau = " << "  " << delta_tau << endl;
	cout << "RPY_END_DESIRED:" << endl;
	cout << RPY_END_DESIRED[0] << "  " << RPY_END_DESIRED[1] << "  " << RPY_END_DESIRED[2] << endl;
	cout << "Pe_DESIRED:" << endl;
	cout << Pe_DESIRED[0] << "  " << Pe_DESIRED[1] << "  " << Pe_DESIRED[2] << endl;

	cout << "Iter_Max = " << "  " << Iter_Max << endl;
	cout << "num_calc = " << "  " << num_calc << endl;
	cout << "particleCount = " << "  " << particleCount << endl;
	cout << "inertGuideCoe = " << "  " << inertGuideCoe << endl;
	cout << "globalGuideCoe = " << "  " << globalGuideCoe << endl;
	cout << "localGuideCoe = " << "  " << localGuideCoe << endl;
	cout << "disturbanceRate = " << "  " << disturbanceRate << endl;
	cout << "disturbanceVelocityCoe = " << "  " << disturbanceVelocityCoe << endl;

	double* result_min = new double[N];
	double* result_max = new double[N];
	double* para_low_bound = new double[N];
	double* para_up_bound = new double[N];

// *******************************************************************************************************************************************************************
	calc_para_range(result_min, result_max);

	for (int i = 0; i<N; i++){
		para_low_bound[i] = result_min[i];
		para_up_bound[i] = result_max[i];
	}
	cout << "result_min:" << endl;
	for(int i=0;i<N;i++){
		cout << result_min[i] << "  ";
	}
	cout << endl;
	cout << "result_max:" << endl;
	for(int i=0;i<N;i++){
		cout << result_max[i] << "  ";
	}
	cout << endl;
// *******************************************************************************************************************************************************************

	//构建算法
	PSO_Algorithm pso(calc_fitness_pso, para_low_bound, para_up_bound, dimension, particleCount, maxSpeed, minSpeed, inertGuideCoe, globalGuideCoe, localGuideCoe);

// *****************************************************************************************************************************************************************
	cout << "算法构建完成，进入迭代搜索过程......" << endl;

	for (int kk = 0; kk < num_calc; kk++){	
		cout << "第" << kk << "次搜索......" << endl;
		start = clock();

        pso.findMin(Iter_Max, bestParticle, disturbanceRate, disturbanceVelocityCoe);

        cout << "best fitness found is: " << "  " << bestParticle._fitness << endl;
		for (int j = 0; j < dimension; j++){
			cout << bestParticle._position[j] << "  ";
		}
		cout << endl;

	    end = clock();   //结束时间
	    cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;
    }

	delete[] para_low_bound;
	delete[] para_up_bound;
	delete[] result_min;
	delete[] result_max;
	delete[] minSpeed;
	delete[] maxSpeed;

}


// cost function 考虑末端位姿, and locus, and base attitude
double calc_fitness_pso(PSO_Particle& particle){
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
	// double* collision_times = new double;

	for (int i = 0; i < N; i++){
		para[i] = particle._position[i];
	}

    double K_a = 1 / 0.0002;  // 1 degree error tolerance
	double K_p = 1 / 0.005;  //  0.002 meter error tolerance for end-effector
    double K_s = 1 / 0.1;  // locus tolerance 0.1 meter.
    double K_b = 1 / 0.0008;  // attitude error tolerance for base is 5 degree
    double K_M = 0;
    double K_t = 0;  //  max allowed motion time is set to 100 seconds, tolerance is 10 seconds.

	double cost_func = 0;

	double delta_xi_end_mod_temp = 0;
	double delta_Pe_end_mod_temp = 0;
	double delta_xi_base_mod_temp = 0;

	forward_kin_2(para, eta_end,  xi_end, Pe, eta_b, xi_b, p_e_initial, locus, delta_xi_b_distrb_max, manipl, T_min);

// ********************************************************  straight_line_locus  ************************************************************************************
    double delta_p_e = 0;
    double straight_line_locus = 0;
    for(int i=0;i<3;i++){
        delta_p_e = p_e_initial[i] - Pe[i];
        straight_line_locus += delta_p_e * delta_p_e;
    }
    straight_line_locus = sqrt(straight_line_locus);
// **********************************************************************************************************************************************************************

	zyx2quaternion(RPY_END_DESIRED[2], RPY_END_DESIRED[1], RPY_END_DESIRED[0], quaternion_end_desired);
	*eta_end_desired = quaternion_end_desired[0];
	for (int i = 0; i<3; i++){
		xi_end_desired[i] = quaternion_end_desired[i + 1];
	}
	delta_var(Pe_DESIRED, eta_end_desired, xi_end_desired, Pe, *eta_end, xi_end, *eta_b, xi_b, delta_eta_end, delta_xi_end, delta_Pe_end, delta_eta_base, delta_xi_base);
	for (int i = 0; i < 3; i++){
		delta_xi_end_mod_temp += delta_xi_end[i] * delta_xi_end[i];
		delta_Pe_end_mod_temp += delta_Pe_end[i] * delta_Pe_end[i];
		delta_xi_base_mod_temp += delta_xi_base[i] * delta_xi_base[i];
	}
	delta_xi_end_mod_temp = sqrt(delta_xi_end_mod_temp); 
    delta_Pe_end_mod_temp = sqrt(delta_Pe_end_mod_temp); 
    delta_xi_base_mod_temp = sqrt(delta_xi_base_mod_temp);

	// ****************************************  RPY error of end-effector,  Pe error  ********************************************************************************

    cost_func = K_a * delta_xi_end_mod_temp 
                + K_p * delta_Pe_end_mod_temp 
                + K_b * (*delta_xi_b_distrb_max) 
                + K_s * fabs((*locus) - straight_line_locus)  
                + K_M * (1 / (*manipl)) 
                + K_t * (*T_min);
				// + (*collision_times);


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
	// delete collision_times;

	return cost_func;

}




// double calc_fitness_tlbo(Student& student){
// 	double* para = new double[N];
// 	double* eta_end = new double;
// 	double* xi_end = new double[3];
// 	double* Pe = new double[3];
// 	double* eta_b = new double;
// 	double* xi_b = new double[3];
// 	double* quaternion_end_desired = new double[4];
// 	double* eta_end_desired = new double;
// 	double* xi_end_desired = new double[3];
// 	double* delta_eta_base = new double;
// 	double* delta_eta_end = new double;
// 	double* delta_xi_base = new double[3];
// 	double* delta_xi_end = new double[3];
// 	double* delta_Pe_end = new double[3];
//     double* p_e_initial = new double[3];
//     double* locus = new double;
//     double* delta_xi_b_distrb_max = new double;
//     double* manipl = new double;
//     double* T_min = new double;
// 	// double* collision_times = new double;
// 	for (int i = 0; i < N; i++){
// 		para[i] = student._position[i];
// 	}
// 	// double K_a = 1 / 0.0087;
//     double K_a = 1 / 0.0002;  // 1 degree error tolerance for end-effector
// 	// double K_p = 1 / 0.005;  //  0.005 meter error tolerance for end-effector
// 	double K_p = 1 / 0.002;  //  0.002 meter error tolerance for end-effector
//     // double K_s = 1 / 0.1;  // locus tolerance 0.1 meter for end-effector.
//     double K_s = 1 / 0.1;  // locus tolerance 0.1 meter for end-effector.
//     double K_b = 1 / 0.0008;  // attitude error tolerance for base is 5 degree
//     double K_M = 0;
//     double K_t = 0;  //  max allowed motion time is set to 100 seconds, tolerance is 10 seconds.
// 	double cost_func = 0;
// 	double delta_xi_end_mod_temp = 0;
// 	double delta_Pe_end_mod_temp = 0;
// 	double delta_xi_base_mod_temp = 0;
// 	forward_kin_2( para,  eta_end,  xi_end,  Pe,  eta_b,  xi_b,  p_e_initial,  locus,  delta_xi_b_distrb_max,  manipl, T_min);
// // ********************************************************  straight_line_locus  ********************************************************************
//     double delta_p_e = 0;
//     double straight_line_locus = 0;
//     for(int i=0;i<3;i++){
//         delta_p_e = p_e_initial[i] - Pe[i];
//         straight_line_locus += delta_p_e * delta_p_e;
//     }
//     straight_line_locus = sqrt(straight_line_locus);
// // ***************************************************************************************************************************************************
// 	zyx2quaternion(RPY_END_DESIRED[2], RPY_END_DESIRED[1], RPY_END_DESIRED[0], quaternion_end_desired);
// 	*eta_end_desired = quaternion_end_desired[0];
// 	for (int i = 0; i<3; i++){
// 		xi_end_desired[i] = quaternion_end_desired[i + 1];
// 	}
// 	delta_var(Pe_DESIRED, eta_end_desired, xi_end_desired, Pe, *eta_end, xi_end, *eta_b, xi_b, delta_eta_end, delta_xi_end, delta_Pe_end, delta_eta_base, delta_xi_base);
// 	for (int i = 0; i < 3; i++){
// 		delta_xi_end_mod_temp += delta_xi_end[i] * delta_xi_end[i];
// 		delta_Pe_end_mod_temp += delta_Pe_end[i] * delta_Pe_end[i];
// 		delta_xi_base_mod_temp += delta_xi_base[i] * delta_xi_base[i];
// 	}
// 	delta_xi_end_mod_temp = sqrt(delta_xi_end_mod_temp);
// 	delta_Pe_end_mod_temp = sqrt(delta_Pe_end_mod_temp);
// 	delta_xi_base_mod_temp = sqrt(delta_xi_base_mod_temp);
//     cost_func = K_a * delta_xi_end_mod_temp 
//                 + K_p * delta_Pe_end_mod_temp 
//                 + K_b * (*delta_xi_b_distrb_max) 
//                 + K_s * fabs(*locus - straight_line_locus)  
//                 + K_M * (1 / (*manipl)) 
//                 + K_t * (*T_min);
// 				// + (*collision_times);
//     delete[] para ;
//     delete eta_end ;
//     delete[] xi_end ;
//     delete[] Pe ;
//     delete eta_b ;
// 	delete[] xi_b ;
// 	delete[] quaternion_end_desired; 
// 	delete eta_end_desired ;
// 	delete[] xi_end_desired ;
// 	delete delta_eta_base ;
// 	delete delta_eta_end ;
// 	delete[] delta_xi_base ;
// 	delete[] delta_xi_end ;
// 	delete[] delta_Pe_end ;
// 	delete[] p_e_initial;
// 	delete locus ;
// 	delete delta_xi_b_distrb_max ;
// 	delete manipl;
// 	delete T_min ;
// 	// delete collision_times;
// 	return cost_func;
// }



// 数组排序，从大到小排序，并返回排序后的数组对应原数组的下标
// void sort_max2min_main(double* a, int length, int* b)
// {
// 	int temp_index;
// 	double temp;
// 	for (int j = 0; j < length; j++)
// 	{
// 		for (int i = 0; i < length - 1 - j; i++)
// 		{
// 			if (a[i] < a[i + 1])
// 			{
// 				temp = a[i];
// 				a[i] = a[i + 1];
// 				a[i + 1] = temp;
// 				temp_index = b[i];
// 				b[i] = b[i + 1];
// 				b[i + 1] = temp_index;
// 			}
// 		}
// 	}
// }














