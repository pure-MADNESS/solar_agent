/*

  ____        _            _____ _  _______    ____ _               
 / ___|  ___ | | __ _ _ __| ____| |/ /  ___|  / ___| | __ _ ___ ___ 
 \___ \ / _ \| |/ _` | '__|  _| | ' /| |_    | |   | |/ _` / __/ __|
  ___) | (_) | | (_| | |  | |___| . \|  _|   | |___| | (_| \__ \__ \
 |____/ \___/|_|\__,_|_|  |_____|_|\_\_|      \____|_|\__,_|___/___/
                                                                    

*/

#ifndef __SOLAREKF_H__
#define __SOLAREKF_H__

#include "ekf.hpp"

using namespace Eigen;

class SolarEKF : public EKF {
    public:
        /**
         * @brief State: [G, P_max]^T
         *        Measurement: [I_dc] measured current at panel terminals
         */
        SolarEKF(double area, double efficiency);

        void set_inputs(double g_api, double v_now);

        /**
         * Prediction:
         * P = A eta G
         */
        Eigen::VectorXd f(const Eigen::VectorXd& x, double dt) override;

        Eigen::MatrixXd F(const Eigen::VectorXd& x, double dt) override;

        Eigen::VectorXd h(const Eigen::VectorXd& x_pred) override;

        Eigen::MatrixXd H(const Eigen::VectorXd& x) override;

    private:
        double _area;       // [m^2]
        double _eff;        // efficiency
        double _g_forecast; // [W/m^2]
        double _v_actual;   // [V] measured DC voltage
};


#endif
