#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


// sets weights for the UKF.
void UKF::SetWeights(){
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i< 2 * n_aug_ + 1; i++) {  //2n+1 weights
      double weight = 0.5/(n_aug_ + lambda_);
      weights_(i) = weight;
    }
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // state vector dimension
  n_x_ = 5;

  // augmented state dimension
  n_aug_ = 7;

  // lidar sensor measurement vector dimension
  n_z_lid_ = 2;

  // radar sensor measurement vector dimension
  n_z_rad_ = 3;

  //define spreading parameter
  lambda_ = 3 - n_x_;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 0.15,  0,  0, 0, 0,
           0, 0.15, 0, 0, 0,
           0,  0, 0.5, 0, 0,
           0,  0, 0, 0.5, 0,
            0, 0, 0, 0, 0.5;

  //P_ = MatrixXd(n_x_, n_x_);


  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Set Weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  SetWeights();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.2;

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

  // Lidar measurement noise covaraince matrix
  R_lid_ = MatrixXd::Zero(2, 2);
  VectorXd r_var_vec_lid = VectorXd(2);
  r_var_vec_lid << (std_laspx_ * std_laspx_), (std_laspy_ * std_laspy_);
  // covariances set to 0;
  R_lid_.diagonal() = r_var_vec_lid;

  // Radar measurement noise covaraince matrix
  R_rad_ = MatrixXd::Zero(n_z_rad_, n_z_rad_);
  VectorXd r_std_vec_rad = VectorXd(n_z_rad_);
  r_std_vec_rad << std_radr_, std_radphi_, std_radrd_;
  r_std_vec_rad = r_std_vec_rad.array().square();
  // covariances set to 0;
  R_rad_.diagonal() = r_std_vec_rad;
}//end UKF constructor


  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

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
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // Initialize state with first measurement
    std::cout << "EKF Initialization... " << std::endl;
    if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];
      float v = 1;
      float yaw = 1;
      float yawrate = 0.1;

      x_ << px, py, v, yaw, yawrate;

    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float px = rho * std::cos(phi);
      float py = rho * std::sin(phi);
      float v = 1;
      float yaw = 1;
      float yawrate = 0.1;

      x_ << px, py, v, yaw, yawrate;
    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    std::cout << "UKF Initialized!" << std::endl;
    return;
  } //if


  //compute the time elapsed between the current and previous measurements
  //delta_t - expressed in seconds
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  /* ### PREDICT State ###*/
  // call the prediction routine
  Prediction(delta_t);

  /* ###  UPDATE State ###*/
  // depending on sensor type, call relevant update routine
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    // update using LIDAR sensor measurement
    UpdateLidar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    // update using Radar sensor measurement
    UpdateRadar(meas_package);
  }

}//end ProcessMeasurement

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

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Generate sigma points
   //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v / yawd * ( sin (yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v / yawd * ( cos(yaw) - cos(yaw + yawd * delta_t) );
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    float angle = x_diff(3);
    x_diff(3) = atan2(sin(angle), cos(angle));

    // while (x_diff(3)> M_PI) x_diff(3) -= 2. * M_PI;
    // while (x_diff(3)<-M_PI) x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

} // end Prediction

void UKF::PredictLidarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S){

  //transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred_.cols(); i++) {  //2n+1 sigma points

      // extract values for better readibility
      double p_x = Xsig_pred_(0,i);
      double p_y = Xsig_pred_(1,i);

      // measurement model to transform sigma points to measurment space
      Zsig(0,i) = p_x; //x
      Zsig(1,i) = p_y; //y
  }

  // calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;

      S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix

  S = S + R_lid_;

} // end PredictLidarMeasurement


void UKF::PredictRadarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S){

  //transform sigma points into measurement space
  for (int i=0; i<Xsig_pred_.cols(); i++){
     float px = Xsig_pred_.col(i)[0];
     float py = Xsig_pred_.col(i)[1];
     float v = Xsig_pred_.col(i)[2];
     float psi = Xsig_pred_.col(i)[3];
     float psi_dot = Xsig_pred_.col(i)[4];

     float px_py_norm = sqrt(px*px + py*py);

     // measurement model to transform sigma points to measurment space
     Zsig.col(i)[0] = px_py_norm; //rho
     Zsig.col(i)[1] = atan2(py, px); //phi
     Zsig.col(i)[2] = (px * cos(psi) * v + py * sin(psi) * v) / px_py_norm; //rho_dot
  }

  // calculate mean predicted measurement
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
     z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  // calculate measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++){
    // predicted measurement difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    float angle = z_diff(1);
    z_diff(1) = atan2(sin(angle), cos(angle));

    S = S + weights_(i) * z_diff * z_diff.transpose() ;
  }
  S = S + R_rad_;
} // end PredictLidarMeasurement


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
  // Initialize sigma points
  // Augmented sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_lid_, 2 * n_aug_ + 1);
  // predicted measurement vector
  VectorXd z_pred = VectorXd(n_z_lid_);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_lid_, n_z_lid_);

  PredictLidarMeasurement(Zsig, z_pred, S);

  // vector of sensor measurements
  VectorXd z = VectorXd(n_z_lid_);
  z = meas_package.raw_measurements_;

  // update the state of the UKF:
  UpdateState(Zsig, z_pred, S, z, meas_package.sensor_type_, NIS_lidar_);
} // End UpdateLidar

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
  // Initialize sigma points
  // Augmented sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_rad_, 2 * n_aug_ + 1);
  // predicted measurement vector
  VectorXd z_pred = VectorXd(n_z_rad_);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_rad_, n_z_rad_);

  PredictRadarMeasurement(Zsig, z_pred, S);

  // vector of sensor measurements
  VectorXd z = VectorXd(n_z_rad_);
  z = meas_package.raw_measurements_;

  // update the state of the UKF:
  UpdateState(Zsig, z_pred, S, z, meas_package.sensor_type_, NIS_radar_);

} // end UpdateRadar


void UKF::UpdateState(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S,
                      VectorXd& z, MeasurementPackage::SensorType& sensor, double& NIS){
  /*
    Updates the state of the UKF using the measurement vector z.
  */

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, z.size());


  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization for radar
    if (sensor == MeasurementPackage::RADAR){
      double angle = z_diff(1);
      z_diff(1) = atan2(sin(angle), cos(angle));
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    double angle = x_diff(3);
    x_diff(3) = atan2(sin(angle), cos(angle));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // S inverse
  MatrixXd S_inv = S.inverse();
  // Kalman gain K;
  MatrixXd K = Tc * S_inv;

  // measurement residual
  VectorXd z_diff = z - z_pred;

  //angle normalization for Radar
  if (sensor == MeasurementPackage::RADAR){
    double angle = z_diff(1);
    z_diff(1) = atan2(sin(angle),cos(angle));
  }

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // calculate NIS
  NIS = z_diff.transpose() * S_inv * z_diff;
  //std::cout << "NIS: " << NIS << std::endl;
} // end UpdateState
