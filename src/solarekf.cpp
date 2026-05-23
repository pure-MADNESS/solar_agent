/*

  ____        _            _____ _  _______    ____ _               
 / ___|  ___ | | __ _ _ __| ____| |/ /  ___|  / ___| | __ _ ___ ___ 
 \___ \ / _ \| |/ _` | '__|  _| | ' /| |_    | |   | |/ _` / __/ __|
  ___) | (_) | | (_| | |  | |___| . \|  _|   | |___| | (_| \__ \__ \
 |____/ \___/|_|\__,_|_|  |_____|_|\_\_|      \____|_|\__,_|___/___/
                                                                    

*/

#include "solarekf.hpp"

SolarEKF::SolarEKF(double area, double efficiency) 
  : EKF(2, 2), _area(area), _eff(efficiency) {

  Q.resize(2, 2);
  R.resize(2, 2);

  Q << 0.1, 0, 
       0,   0.1;

  // current sensor uncertainty
  R << 0.1, 0,
        0, 0.01; 
}

void SolarEKF::set_inputs(double v_now, double p_actual) {
  _v_actual = std::max(v_now, 0.1);
  _p_actual = p_actual;
}

VectorXd SolarEKF::f(const VectorXd& x, double dt) {
    VectorXd x_new(2);
    
    x_new(0) = x(0); 
    double p_irr = _area * _eff * x(0);
    x_new(1) = 0.95 * x(1) + 0.05 * p_irr;

    return x_new;
}

MatrixXd SolarEKF::F(const VectorXd& x, double dt) {
    MatrixXd Fj(2, 2);
    
    Fj << 1.0,      0.0,
          0.05 * _area * _eff,      0.95;

    return Fj;
}
VectorXd SolarEKF::h(const VectorXd& x_pred) {
  VectorXd z_pred(2);

  z_pred(0) = x_pred(0);
  double p_expected = std::min(_p_actual, x_pred(1));
  z_pred(1) = p_expected / _v_actual;
  
  return z_pred;
}

MatrixXd SolarEKF::H(const VectorXd& x) {
  MatrixXd Hj(2, 2);
  
  double dI_dP = 0.0;
  // Se siamo in saturazione (la richiesta supera o eguaglia la capacità stimata),
  // allora una variazione della capacità massima altera la corrente misurata!
  if (_p_actual >= x(1)) {
    dI_dP = 1.0 / _v_actual;
  } else {
    dI_dP = 0.0; // Sotto il limite di saturazione, la corrente dipende solo dal carico
  }

  Hj << 1.0, 0.0,
        0.0, dI_dP;
  
  return Hj;
}