/*

  ____        _            _____ _  _______    ____ _               
 / ___|  ___ | | __ _ _ __| ____| |/ /  ___|  / ___| | __ _ ___ ___ 
 \___ \ / _ \| |/ _` | '__|  _| | ' /| |_    | |   | |/ _` / __/ __|
  ___) | (_) | | (_| | |  | |___| . \|  _|   | |___| | (_| \__ \__ \
 |____/ \___/|_|\__,_|_|  |_____|_|\_\_|      \____|_|\__,_|___/___/
                                                                    

*/

#include "solarekf.hpp"

SolarEKF::SolarEKF(double area, double efficiency) 
  : EKF(2, 1), _area(area), _eff(efficiency) {

  Q.resize(2, 2);
  R.resize(1, 1);

  Q << 0.01, 0, 
       0,   1.0;

  // current sensor uncertainty
  R << 0.01; 
}

void SolarEKF::set_inputs(double g_api, double v_now) {
  _g_forecast = g_api;
  _v_actual = std::max(v_now, 0.1);
}

VectorXd SolarEKF::f(const VectorXd& x, double dt) {
    VectorXd x_new(2);
    double alpha = 0.98; 
    
    x_new(0) = alpha * x(0) + (1.0 - alpha) * _g_forecast;  
    x_new(1) = _area * _eff * x_new(0); 

    return x_new;
}

MatrixXd SolarEKF::F(const VectorXd& x, double dt) {
    MatrixXd Fj(2, 2);
    
    double alpha = 0.98;
    double k = _area * _eff;

    Fj << alpha, 0.0,
          k,     0.0; 
    
    return Fj;
}
VectorXd SolarEKF::h(const VectorXd& x_pred) {
  VectorXd z_pred(1);
  
  double p_max = x_pred(1);
  z_pred(0) = p_max / _v_actual;
  
  return z_pred;
}

MatrixXd SolarEKF::H(const VectorXd& x) {
  MatrixXd Hj(1, 2);
  
  // dh/dG = 0
  // dh/dP_max = 1 / V_actual
  Hj << 0, (1.0 / _v_actual);
  
  return Hj;
}