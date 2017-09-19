#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0,      0,
              0,    0.0009, 0,
              0,    0,      0.09;

  //LIDAR H matrix
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ << 1, 1, 0, 0,
         1, 1, 0, 0, 
         1, 1, 1, 1;

  noise_ax = 9; //given value
  noise_ay = 9; //given value

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
  {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 5.2, 0;  //vx and vy are important for RSME

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Use first measurements to initialize the state vectors and covariance
      // matrices.
      // set ekf_.x_(0) to ro*cos(theta)  // convert polar to Cartisan coords
      // set ekf_.x_(1) to ro*sin(theta)
      ekf_.x_(0) = measurement_pack.raw_measurements_(0)*
                      cos(measurement_pack.raw_measurements_(1));
      ekf_.x_(1) = measurement_pack.raw_measurements_(0)*
                      sin(measurement_pack.raw_measurements_(1));
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      // Use first measurements to initialize the state vectors and covariance
      // matrices.
      // set ekf_.x_(0) to x
      // set ekf_.x_(1) to y
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
    }

    // 4X4 state transition matrix
    ekf_.F_ = MatrixXd(4,4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    // 4X4 state covariance matrix
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;

    previous_timestamp_ = measurement_pack.timestamp_;  //VERY IMPORTANT

    // print the output
    cout << "x_init = " << endl << ekf_.x_ << endl;
    cout << "P_init = " << endl << ekf_.P_ << endl << endl;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Lesson 5 section 8
  // Need to compute delta t
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //secs
  previous_timestamp_ = measurement_pack.timestamp_;  

  // Set up some easy conventions
  float dt_2 = dt * dt;
  float dt_3 = dt * dt * dt;
  float dt_4 = dt * dt * dt * dt;

  // Modify the F matrix with the updated time delta (Lesson 5 section 8)
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  // Set the process covariance matrix Q (Lesson 5 section 9)
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
             dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_; 
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << endl << ekf_.x_ << endl;
  cout << "P_ = " << endl << ekf_.P_ << endl << endl;
}
