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

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check validity of inputs
  if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
    std::cout << "Invalid estimation or ground truth data!" << std::endl;
    return rmse;
  }
  // accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];

    // square the residuals, component-wise (Hadamard) product
    residual = residual.array() * residual.array();

    rmse += residual;
  }

  // calculate mean
  rmse = rmse / estimations.size();

  // square root
  rmse = rmse .array().sqrt();

  return rmse;

}
