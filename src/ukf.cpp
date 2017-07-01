#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
#define EPS 0.001

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
  std_a_ = 10; //30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 10; //30

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
  n_aug_=7; //Add 2 more augumented state variable x_ to accomodate process noise
  lambda_= 3-n_aug_; //Empirical Formula
  n_sigpts_ = 2*n_aug_+1;
  is_initialized_ = false;

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

  //Initialize weights for sigma points
  //set weights
 weights_ = VectorXd::Zero(n_sigpts_);
  for (int i=1; i < n_sigpts_; i++)
    weights_(i) = 0.5/(lambda_+n_aug_);
  weights_(0)= lambda_/(lambda_+n_aug_);

 //Initialize Predicted Sigma Points
 Xsig_pred_ = MatrixXd::Zero(n_x_,n_sigpts_);
  
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
  double dt = (meas_package.timestamp_- prev_timestamp_us_)/1e6;
  //Update previous timestamp, forgetting this leads nan values
  prev_timestamp_us_ = meas_package.timestamp_;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "Initializing Matrices: " << endl;

     //Initialize P,Q,R,F matrices
     // initial state vector to have 5 states (px,py,v,psi,psi_dot)
      x_ = VectorXd::Zero(n_x_);
      
      // initial covariance matrix to covariance for 5 states
      P_ = MatrixXd::Identity(n_x_, n_x_);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float r0 = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float rdot = meas_package.raw_measurements_[2];
      x_ << r0*cos(theta), r0*sin(theta),rdot,0,0 ; //* Need to figure our velocity
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];
      x_ << px,py,0,0,0;
    }

     // Deal with the special case initialisation problems
    if (fabs(x_(0)) < EPS and fabs(x_(1)) < EPS){
		x_(0) = EPS;
		x_(1) = EPS;
    }
  
    is_initialized_ = true;
    cout << "End of Initialization";
    cout << "Initialized X state" << x_ << endl;

    return; //No prediction for first intialization
  } //End of Init

    Prediction(dt);
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      UpdateLidar(meas_package);
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
  Tools tools;
 /*******************************************************************************
 * Calculate Augumented Sigma Points
 ******************************************************************************/
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigpts_);
  
  //init augmented mean state
  x_aug.head(n_x_) = x_;
  //cout << "Augumented Mean State " << x_aug << endl;

  //init augmented covariance matrix
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create and compute square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //init augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug  + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug  - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
 /*******************************************************************************
 * Compute Predicted Sigma Points
 ******************************************************************************/
  float px,py,v,psi,psi_dot,vak,vpsik,dt,dt2;
  VectorXd pred_x(n_x_),ext_x(n_x_);
  
  for (int i=0; i < n_sigpts_; i++)
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
    
      if (fabs(psi_dot) > EPS){
        pred_x << v/psi_dot*(sin(psi+psi_dot*dt) - sin(psi)),
        v/psi_dot*(-cos(psi+psi_dot*dt) +cos(psi)),0,
        psi_dot*dt,0;
        }
        else {
        pred_x << v*cos(psi)*dt,v*sin(psi)*dt,0,psi_dot*dt,0;
        }
      //cout << "i= " << i << endl;
      //cout << pred_x << endl;

      ext_x << 0.5*dt2*cos(psi)*vak,0.5*dt2*sin(psi)*vak, dt*vak,0.5*dt2*vpsik,dt*vpsik;
      pred_x += ext_x;
      pred_x += Xsig_aug.col(i).head(n_x_); //Update  current state
      //cout << pred_x.size() << endl;
    //Write into output Matrix
     Xsig_pred_.col(i) << pred_x;

  }

 /*******************************************************************************
 * Compute Predicted State's Mean and Covariance 
 ******************************************************************************/
  //predict state mean, sum of weights =1
  x_.fill(0.0);
  for (int i=0; i < n_sigpts_; i++)
    x_ += weights_(i)*Xsig_pred_.col(i);

  //predict state covariance matrix
  P_.fill(0);
  VectorXd x_diff(n_x_);
  for (int i=0; i < n_sigpts_; i++)
  {
    //cout << "x" << endl << x_ << endl;
    x_diff = (Xsig_pred_.col(i)-x_);
    x_diff(3) = tools.NormalizedAngle(x_diff(3));
    //cout << "x_diff" << endl << x_diff << endl;
    P_ += weights_(i)*x_diff*x_diff.transpose();
  }

  //cout << "P matrix: " << endl << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  Also calculate the lidar NIS.
  */
 /*******************************************************************************
 * Transform sigma points into measurement space
 ******************************************************************************/
  MatrixXd Xsig_pred(n_aug_,n_sigpts_);
  Xsig_pred = Xsig_pred_;
  
  int n_z =2; //Only 3 measurements for LIDAR sensor
  MatrixXd Zsig = MatrixXd::Zero(n_z,n_sigpts_);

  //transform sigma points into LIDAR measurement space
  for (int i = 0; i < n_sigpts_; i++) { 

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model (px,py for LIDAR)
    Zsig(0,i) = p_x;                       
    Zsig(1,i) = p_y; 
  }

 /*******************************************************************************
 * Update state variables using Measurement Data by computing S and Tc matrices(-> K ->x)
 ******************************************************************************/
  VectorXd z = meas_package.raw_measurements_;
  compute_SandT_UpdateState(Zsig,z,meas_package);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  Also calculate the radar NIS.
  */
 /*******************************************************************************
 * Transform sigma points into measurement space
 ******************************************************************************/
  MatrixXd Xsig_pred(n_aug_,n_sigpts_);
  Xsig_pred = Xsig_pred_;
  
  int n_z =3; //Only 3 measurements for RADAR sensor
  MatrixXd Zsig(n_z,n_sigpts_);

  //transform sigma points into RADAR measurement space
  for (int i = 0; i < n_sigpts_; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // RADAR measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }


 /*******************************************************************************
 * Update state variables using Measurement Data by computing S and Tc matrices(-> K ->x)
 ******************************************************************************/
  VectorXd z = meas_package.raw_measurements_;
  compute_SandT_UpdateState(Zsig,z,meas_package);

}

void UKF::compute_SandT_UpdateState(MatrixXd Zsig, VectorXd z,MeasurementPackage meas_package)
{
  //Inputs: Zsig, z
  //input class variables: x_
  //Updated class variables: x_, P_
 
   int n_z =z.size(); //Only 3 measurements for RADAR sensor
   Tools tools;

  //Mean predicted measurement from Model(required for S computation and Kalman Update)
  VectorXd zmean_pred = VectorXd(n_z);
  zmean_pred.fill(0.0);
  for (int i=0; i < n_sigpts_; i++) {
      zmean_pred = zmean_pred + weights_(i) * Zsig.col(i);
  }

  //cout << "Current Mean Predicted State" << endl << zmean_pred << endl;
  //Get current state
  VectorXd x = x_;

 /*******************************************************************************
 * Compute S Matrix
 ******************************************************************************/
  //Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  
  for (int i = 0; i < n_sigpts_; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - zmean_pred;

    //angle normalization
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      z_diff(1) = tools.NormalizedAngle(z_diff(1)); //Normalize phi of RADAR diff
    }

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //cout << "S Matrix: " << endl << S << endl;

  //add measurement noise covariance matrix
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    S = S + R_radar_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    S = S + R_lidar_;
  }

 /*******************************************************************************
 * Compute Tc Matrix
 ******************************************************************************/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - zmean_pred;
    
    //angle normalization
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      z_diff(1) = tools.NormalizedAngle(z_diff(1)); //Normalize phi of RADAR diff
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x; //main compute
    
    //angle normalization
    x_diff(3) =   tools.NormalizedAngle(x_diff(3)); //Normalize psi

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  //cout << "Tc Matrix: " << endl << Tc << endl;

  /******************************************************************************
  * Update State using S and T matrices
  *******************************************************************************/
  VectorXd zmean_diff = z - zmean_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
    // Angle normalization
    zmean_diff(1) =   tools.NormalizedAngle(zmean_diff(1));
  }
  UpdateState_from_SnT(S,Tc,zmean_diff);
  cout << "Updated State from Measurement" << endl;
}

void UKF::UpdateState_from_SnT(MatrixXd S, MatrixXd Tc, VectorXd zmean_diff)
{
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  x_ = x_ + K * zmean_diff;
  P_ = P_ - K*S*K.transpose();

  cout << "x ="  << endl << x_ << endl;
  cout << "P = " << endl << P_ << endl;

}