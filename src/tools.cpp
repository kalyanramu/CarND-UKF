#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  int x_size = estimations[0].size();
  VectorXd rmse = VectorXd::Zero(x_size);
  for(int i=0; i < estimations.size(); ++i){
        // ... your code here
  		VectorXd residual = estimations[i] - ground_truth[i];
  		residual = residual.array()*residual.array();
  		rmse += residual;
	}

	rmse = rmse/estimations.size(); //Divide by number of measurements
	rmse = rmse.array().sqrt();
	return rmse;
}

double Tools::NormalizedAngle(double Angle)
{
     //angle normalization
    // while (Angle > M_PI) Angle -=2.*M_PI;
    // while (Angle <-M_PI) Angle +=2.*M_PI;
    // return Angle;

  float M_2PI = 2*M_PI;
  Angle = fmod(Angle + M_PI,M_2PI);
  if (Angle < 0)
        Angle += M_2PI;
  return Angle - M_PI;
}