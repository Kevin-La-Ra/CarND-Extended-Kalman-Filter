#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  VectorXd rmse(4);

  rmse << 0,0,0,0;

  // check the input validity
  if (estimations.size() == 0 || estimations.size() != ground_truth.size())
  {
    cout << "Invalid estimation size or ground_truth data" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i)
  {
    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{  MatrixXd Hj(3,4);  //recover state parameters  float px = x_state(0);  float py = x_state(1);  float vx = x_state(2);  float vy = x_state(3);  float t1 = px*px + py*py;  float t2 = sqrt(px*px + py*py);  float t3 = (t1*t2);  //check for division by zero  if(fabs(t1) < 0.0001)
  {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;    return Hj;  }  else
  {    //compute the Jacobian matrix    Hj << (px/t2), (py/t2), 0, 0, -(py/t1), (px/t1), 0, 0,         py*(vx*py - vy*px)/t3, px*(px*vy - py*vx)/t3, px/t2, py/t2;  }
  return Hj;}