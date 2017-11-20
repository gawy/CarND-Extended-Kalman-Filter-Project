#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  cout << "KalmanFilter::Predict" << std::endl;
  /**
    * predict the state
  */
  x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;

  cout << "Predicted x: " << x_ << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  cout << "KalmanFilter::Update" << std::endl;
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K =  P_ * Ht * Si;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

	//new state
	x_ = x_ + (K * y);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  cout << "KalmanFilter::UpdateEKF" << std::endl;
  /**
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd h_x(3);
  float sqrp = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double phi = atan2(x_(1), x_(0));
  h_x << sqrp,
         phi,
         (x_(0)*x_(2) + x_(1)*x_(3)) / sqrp;

  // cout << "h_x: " << h_x << std::endl;

  VectorXd y = z - h_x;


  if (y(1) < -M_PI || y(1) > M_PI) {
    float sign = y(1) > 0 ? -1 : 1;
    int whole = abs((y(1) - sign*M_PI) / (2*M_PI));
    float new_phi = y(1) + sign*whole*2*M_PI;
    cout << "phi in y: " << y(1) << ", adjusted to: "<< new_phi << std::endl;
    y(1) = new_phi;
  }

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  //new state
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}
