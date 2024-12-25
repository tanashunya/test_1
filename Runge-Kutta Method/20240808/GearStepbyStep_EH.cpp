#include <BasicDefines.hpp>

#include "GearStepbyStep_EH.hpp"




/*//=======================================================
  // ●　１ステップ計算
  //=======================================================*/
void GearStepbyStep_EH::compute_step(){
	const double gear_t = 2.0/3.0 * delta_t;
	//const double next_time = time+delta_t;
	double temp_x1 = pos_x1[t_step-1];
	double temp_x2 = pos_x2[t_step-1];
	double temp_v1 = vel_x1[t_step-1];
	double temp_v2 = vel_x2[t_step-1];

	temp_save_results.clear();

	const int solved_size = 2;//4;
	Eigen::MatrixXd matK = Eigen::MatrixXd::Zero(solved_size, solved_size);

	/* [0]:x, [1]:v */
	matK(0, 0) = 1.0;
	matK(0, 1) = -1.0 * gear_t;

	matK(1, 0) = gear_t / mass1 * (spring_k1);
	matK(1, 1) = 1.0 + gear_t / mass1 * (damper_c1);

	/* 右辺ベクトル作成 */
	Eigen::VectorXd vecB;// = Eigen::VectorXd::Zero(solved_size);
	calcRightHand(vecB, time, temp_x1, temp_x2, temp_v1, temp_v2);

	/* 解く */
	Eigen::VectorXd temp_results = matK.inverse()*vecB;

	/*cout << "vib_temp " << time << ", " << t_step << endl <<temp_results << endl;
	cout << "----" << endl;
	cout <<vecB << endl;
	cout << "???" << endl;
	getchar();*/

	/* 一次結果に一式を保存 */
	temp_save_results.clear();
	temp_save_results.resize(10);
	
	double x0 = this->disp_func(time);
	temp_save_results[0] = temp_results(0);/*変位1*/
	temp_save_results[1] = 0;//temp_results(2);/*変位2*/
	temp_save_results[2] = temp_results(1);/*速度1*/
	temp_save_results[3] = 0;//temp_results(3);/*速度2*/
	temp_save_results[4] = temp_results(0)-x0;/*真の変位1*/
	temp_save_results[5] = 0;//temp_results(2)-x0;/*真の変位2*/
}

/*//=======================================================
  // ●　更新結果の確定
  //=======================================================*/
void GearStepbyStep_EH::compute_fix(){
	pos_x1.push_back(temp_save_results[0]);	
	pos_x2.push_back(temp_save_results[1]);
	vel_x1.push_back(temp_save_results[2]);	
	vel_x2.push_back(temp_save_results[3]);
	
	//double x0 = this->disp_func(time);
	true_x1.push_back( temp_save_results[4] );/*真の変位1*/
	true_x2.push_back( temp_save_results[5] );/*真の変位2*/

	time += delta_t;
	t_step++;
}


/*//=======================================================
  // ●　右辺作成
  //=======================================================*/
void GearStepbyStep_EH::calcRightHand(Eigen::VectorXd& vecB, double t, double x1, double x2, double v1, double v2){
	const double gear_t = 2.0/3.0 * delta_t;
	const int solved_size = 2;//4;
	vecB = Eigen::VectorXd::Zero(solved_size);
#endif

	const double input_x = disp_func(t);
	const double input_v = disp_vel_func(t);
	//cout << input_x << ", " <<input_v << endl;

	vecB(0) = 4.0/3.0*pos_x1[t_step-1] - 1.0/3.0*pos_x1[t_step-2];
	vecB(1) = gear_t/mass1 * (damper_c1*input_v + spring_k1*input_x) + 4.0/3.0*vel_x1[t_step-1] - 1.0/3.0*vel_x1[t_step-2];
	/* 電磁力の項 */
	if(mag_f_active){
		vecB(1) += gear_t / mass1 * mag_force;
	}
}



/*//=======================================================
  // ●　変位の時間変化関数
  //=======================================================*/
double GearStepbyStep_EH::disp_func(double tt){
	return vib_amp*sin(wfreq * tt) + noise*sin(wfreq/3.0 * tt);
}

/*//=======================================================
  // ●　変位速度の時間変化関数
  //=======================================================*/
double GearStepbyStep_EH::disp_vel_func(double tt){
	return vib_amp*wfreq*cos(wfreq * tt) + noise*wfreq*cos(wfreq/3.0 * tt) / 3.0;
}

/*//=======================================================
  // ●　変位加速度の時間変化関数
  //=======================================================*/
double GearStepbyStep_EH::disp_accel_func(double tt){
	return -1.0*vib_amp*wfreq*wfreq*sin(wfreq * tt) - noise*wfreq*wfreq*sin(wfreq/3.0 * tt) / 3.0/3.0;
}

