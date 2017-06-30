#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5; //Total number of states in state vector
  n_aug_=2; //Add 2 more augumented state variable to accomodate process noise
  lambda_= 3-n_aug_; //Empirical Formula
  //int n_sig = 2*n_aug+1;

  //Initialize P,Q,R,F matrices
  // initial state vector to have 5 states (px,py,v,psi,psi_dot)
  x_ = VectorXd(5);
  // initial covariance matrix to covariance for 5 states
  P_ = MatrixXd(5, 5);

  //create and initialize R matrices (sensor meas uncertainity matrices) for RADAR and LIDAR 
  //Rradar //rho,phi,pho_dot
  R_radar_ = MatrixXd(3,3);
  R_radar_ <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;

  //Rlidar //px,py
  R_lidar_ = MatrixXd(2,2);
  R_lidar_ <<    std_laspx_*std_laspx_, 0,
               0,std_laspy_*std_laspy_;


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    //UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column


  
 /*******************************************************************************
 * Calculate Augumented Sigma Points
 ******************************************************************************/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);


  //init augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //init augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create and compute square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //init augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
 /*******************************************************************************
 * Compute Predicted Sigma Points
 ******************************************************************************/
  float px,py,v,psi,psi_dot,vak,vpsik,dt,dt2;
  VectorXd pred_x(n_x_),ext_x(n_x_);
  MatrixXd Xsig_pred(n_aug_,2*n_aug_+1); //n_sig = 2*n_aug+1
  
  for (int i=0; i < 2*n_aug_+1; i++)
  {
      //Get state variables
      px = Xsig_aug(0,i);
      py = Xsig_aug(1,i);
      v  = Xsig_aug(2,i);
      psi = Xsig_aug(3,i);
      psi_dot = Xsig_aug(4,i);
      vak = Xsig_aug(5,i);
      vpsik = Xsig_aug(6,i);
      dt = delta_t;
      dt2 = dt*dt;
      
      //Predict current sigma point
    
      if (fabs(psi_dot) > 0.001){
        pred_x << v/psi_dot*(sin(psi+psi_dot*dt) - sin(psi)),
        v/psi_dot*(-cos(psi+psi_dot*dt) +cos(psi)),0,
        psi_dot*dt,0;
        }
        else {
        pred_x << v*cos(psi)*dt,v*sin(psi)*dt,0,psi_dot*dt,0;
        }
        
      ext_x << 0.5*dt2*cos(psi)*vak,0.5*dt2*sin(psi)*vak, dt*vak,0.5*dt2*vpsik,dt*vpsik;
      pred_x += ext_x;
      pred_x += Xsig_aug.col(i).head(n_x_);

    
    //Write into output Matrix
     Xsig_pred.col(i) << pred_x;

  }
  //Assign local variables to class variables
  Xsig_pred_ = Xsig_pred;

 /*******************************************************************************
 * Predict State's Mean and Covariance 
 ******************************************************************************/
  //Init variables
  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);
  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  //set weights
  int w_size = 2*n_aug_+1;
  weights.fill(0);
  for (int i=1; i < w_size; i++)
    weights(i) = 0.5/(lambda_+n_aug_);
  weights(0)= lambda_/(lambda_+n_aug_);

   //predict state mean, sum of weights =1
  x.fill(0.0);
  for (int i=0; i < w_size; i++)
    x += weights(i)*Xsig_pred.col(i);

  //predict state covariance matrix
  P.fill(0);
  VectorXd x_diff(n_x_);
  for (int i=0; i < w_size; i++)
  {
    x_diff = (Xsig_pred.col(i)-x);
    
    //check if angle in state vector is between -pi and pi
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P += weights(i)*x_diff*x_diff.transpose();
    
  }

  //Assign local variables to class variables
  x_ = x;
  P_ = P;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  /*******************************************************************************
 * Transform sigma points into measurement space
 ******************************************************************************/
  MatrixXd Xsig_pred(n_aug_,2*n_aug_+1);
  Xsig_pred = Xsig_pred_;
  
  int n_z =3; //Only 3 measurements for RADAR sensor
  MatrixXd Zsig(n_z,2*n_aug_+1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_radar_;

 /*******************************************************************************
 * Update state variables using Measurement Data (-> K ->x)
 ******************************************************************************/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //Update local variables using class variables
  VectorXd x = x_; 
  MatrixXd P = P_;

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();

  //Update class variables using local variables
  x_ = x;
  P_ = P;
}
