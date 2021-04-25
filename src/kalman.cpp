/*
 * kalman.cpp
 *
 * Author: Iacopo Sprenger
 *
 */

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <math.h>

#include "ros/ros.h"
#include "sensor_msgs/NavSatFix.h"
#include "sensor_msgs/FluidPressure.h"

#include "geometry_msgs/Vector3Stamped.h"

/*
 * CONSTANTS
 */

#define g 		9.81
#define M 		0.2897
#define R		8.3145
#define T		298.15

#define sigma_z_gps	5.0
#define sigma_z_baro	5.0

#define sigma_a		1.0
#define sigma_p0	0.01
#define sigma_k		1e-9


/*
 * MACROS
 */

#define DIAG6(d0, d1, d2, d3, d4, d5)        \
	 (d0),    0,    0,    0,    0,    0, \
	    0, (d1),    0,    0,    0,    0, \
	    0,    0, (d2),    0,    0,    0, \
	    0,    0,    0, (d3),    0,    0, \
	    0,    0,    0,    0, (d4),    0, \
	    0,    0,    0,    0,    0, (d5)

#define DIAG3(d0, d1, d2)   \
	 (d0),    0,    0,  \
	    0, (d1),    0,  \
	    0,    0, (d2)




/*
 * TYPEDEFS
 */

typedef Eigen::Matrix<double, 6, 6> Mat66;

typedef Eigen::Matrix<double, 6, 3> Mat63;

typedef Eigen::Matrix<double, 6, 1> Mat61;

typedef Eigen::Matrix<double, 1, 6> Mat16;

typedef Eigen::Matrix<double, 1, 1> Mat11;

typedef Eigen::Matrix<double, 3, 3> Mat33;


/*
 * GLOBAL VARIABLES
 */

static Mat66 P_tilde;

static Mat61 X_tilde;

static Mat66 P_hat;

static Mat61 X_hat;

static Mat11 R_gps;

static Mat11 R_baro;

static Mat33 Q;

static Mat66 F;

static Mat63 G;

static double last_time;

static Mat66 I;

static ros::Publisher height_pub; 

static bool first;

void initialize(void) {

	X_tilde << 516.0, 0, 0, 9620000, M/(R*T), 516.0;

	P_tilde << DIAG6(25.0, 0.25, 0.25, 25.0, 1e-13, 25.0);
	
	I << DIAG6(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);


	X_hat << X_tilde;

	P_hat << P_tilde;

	R_gps << sigma_z_gps*sigma_z_gps;

	R_baro << sigma_z_baro*sigma_z_baro;

	Q << DIAG3(sigma_a*sigma_a, sigma_p0*sigma_p0, sigma_k*sigma_k);

	F <<     0,  1,  0,  0,  0,  0, 
		 0,  0,  0,  0,  0,  0, 
		 0,  0,  1,  0,  0,  0, 
		 0,  0,  0,  0,  0,  0, 
		 0,  0,  0,  0,  0,  0, 
		 0,  0,  0,  0,  0,  0;

	G <<     0,  0,  0,   
		 0,  0,  0,   
		 1,  0,  0,   
		 0,  1,  0,   
		 0,  0,  1,   
		 0,  0,  0; 

	first = true;


}

void predict(double time) {
	double dt = time - last_time;
	if(dt <= 0) {
		return;
	}

	ROS_INFO("dt: %lf", dt);

	//Create the discrete matrix
	Eigen::Matrix<double, 12, 12> A;
	A.block<6, 6>(0, 0) << -F;
	A.block<6, 6>(6, 0) << DIAG6(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); 
	A.block<6, 6>(6, 6) << F.transpose();
	A.block<6, 6>(0, 6) << G*Q*G.transpose();
	Eigen::Matrix<double, 12, 12> B = A.exp();

	Mat66 PHI;
	PHI << B.block<6, 6>(6, 6).transpose();
	Mat66 Q_w;
	Q_w << PHI*B.block<6, 6>(0, 6);
	
	//Compute step
	X_tilde << PHI*X_hat;
	P_tilde << PHI*P_hat*PHI.transpose() + Q_w;

	last_time = time;

}

void update_baro(double p) {
	double h = X_tilde(0, 0);
	double p0 = X_tilde(3, 0);
	double k = X_tilde(4, 0);
	double h0 = X_tilde(5, 0);

	double p_est = p0*exp(k*g*(h0-h));

	Mat16 H;
	H << -k*g*p_est, 0.0, 0.0, p_est/p0, g*(h0-h)*p_est, k*g*p_est;

	Mat61 K;
	K << P_tilde*H.transpose()*(H*P_tilde*H.transpose() + R_baro).inverse();

	X_hat << X_tilde + K*(p - p_est);
	
	P_hat << (I - K*H)*P_tilde;

	ROS_INFO("baro: %lf", p);


}

void update_gnss(double z) {
	Mat16 H;
	H << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	Mat61 K;
	K << P_tilde*H.transpose()*(H*P_tilde*H.transpose() + R_gps).inverse();

	X_hat << X_tilde + K*(z - H*X_tilde);

	P_hat << (I - K*H)*P_tilde;

	ROS_INFO("gnss: %lf", z);
}

void display(void) {
	ROS_INFO("[%lf] z: %g", last_time, X_tilde(0, 0));
}

void baro_callback(const sensor_msgs::FluidPressure::ConstPtr& data) {
	double pressure = data->fluid_pressure;
	double time = data->header.stamp.nsec/1e9 + data->header.stamp.sec;

	if(first) {
		last_time = time;
		first = false;
		return;
	}

	predict(time);

	update_baro(pressure);

	geometry_msgs::Vector3Stamped pose;
	pose.vector.z = X_tilde(0, 0);
	pose.header.stamp.sec = time;

	height_pub.publish(pose);

	display();

}

void gnss_callback(const sensor_msgs::NavSatFix::ConstPtr& data) {
	double altitude = data->altitude;
	double time = data->header.stamp.nsec/1e9 + data->header.stamp.sec;

	if(first) {
		last_time = time;
		first = false;
		return;
	}


	predict(time);

	update_gnss(altitude);

	geometry_msgs::Vector3Stamped pose;
	pose.vector.z = X_tilde(0, 0);
	pose.header.stamp.sec = time;

	height_pub.publish(pose);

	display();


}

int main(int argc, char ** argv) {
	ros::init(argc, argv, "kalman");

	initialize();

	
	ros::NodeHandle n;

	ros::Subscriber baro_sub = n.subscribe("/iris_1/imu/atm_pressure", 1000, baro_callback);

	ros::Subscriber gnss_sub = n.subscribe("/iris_1/global_position/raw/fix", 1000, gnss_callback);
	
	height_pub = n.advertise<geometry_msgs::Vector3Stamped>("/height", 1000);

	ros::spin();
	
	return 0;
}










