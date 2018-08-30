/*
    USRP_Software_defined_radar is a software for real time sampling, processing, display and storing
    Copyright (C) 2018  Jonas Myhre Christiansen <jonas-myhre.christiansen@ffi.no>
	
    This file is part of USRP_Software_defined_radar.

    USRP_Software_defined_radar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    USRP_Software_defined_radar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with USRP_Software_defined_radar.  If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef TRACKER_H
#define TRACKER_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "detection.h"

using namespace boost::numeric::ublas;

class Tracker
{
private:
    vector<double> x_minus;
    vector<double> x_plus;

    matrix<double> P_minus;
    matrix<double> P_plus;

    double sigmaM, alpha;

    double T, t, t_old;

    matrix<double> F();
    matrix<double> F(double &dt);
    matrix<double> H();
    matrix<double> Q();
    matrix<double> Q(double &dt);
    matrix<double> R(double SNR, double dr, double dv);
    matrix<double> matrixInverse(matrix<double> &m);

public:
    Tracker(double t0, double r0, double v0, double a0, double snr0, double alpha0, double P_0, double P_max, double A_max);

    void motionUpdate(double Tn, bool updateWihtoutDetection);
    void informationUpdate(double rn, double an, double vn, double snrn, double dr, double dv);
    matrix<double> covarianceStateError(double &snrn, double &dr, double &dv, double &dt);
    matrix<double> covarianceMeasureError(double &snrn, double &dr, double &dv);
    void trackPrediction(double &rn, double &an, double &vn, double &snrn);
    void trackEstimate(double &rn, double &an, double &vn, double &snrn);
    vector<double> covariance();
    vector<double> covariancePlus();
    vector<double> distanceInSigmas(double r, double v, double SNR, double dr, double dv);
    double trace(matrix<double>& A);

    double getTimeDiff() {
        return T;
    }

    bool hasBeenUpdated;
    int updatesWithoutDetections;
    bool first_update;

    Detection closest_detection;
};

#endif // TRACKER_H
