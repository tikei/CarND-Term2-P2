#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  // Lidar measurement noise covaraince matrix
  MatrixXd R_lid_;

  // Radar measurement noise covaraince matrix
  MatrixXd R_rad_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Dimension of lidar measurement vector
  int n_z_lid_;

  ///* Dimension of radar measurement vector
  int n_z_rad_;

  ///* Sigma point spreading parameter
  double lambda_;

  // Normalized Innovation Squared (NIS) for LIDAR and RADAR
  double NIS_lidar_;
  double NIS_radar_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */

  void PredictLidarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S);
  /*
  Transforms sigma points into LIDAR measurement space.
  Predicts lidar measurment, before the actual measurement is taken
  into account during the update.
  */

  void PredictRadarMeasurement(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S);
  /*
  Transforms sigma points into RADAR measurement space.
  Predicts radar measurment, before the actual measurement is taken
  into account during the update.
  */

  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  // Update the state the of the UKF
  void UpdateState(MatrixXd& Zsig, VectorXd& z_pred, MatrixXd& S,
                    VectorXd& z, MeasurementPackage::SensorType& sensor, double& NIS);


  private:
    void SetWeights();
};

#endif /* UKF_H */
