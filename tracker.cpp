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
#include "tracker.h"

Tracker::Tracker(double t0, double r0, double v0, double a0, double snr0, double alpha0, double P_0, double P_max, double A_max)
{

    sigmaM = std::sqrt( (std::pow(A_max,2)/3)*(1 + 4*P_max - P_0) );
    alpha = alpha0;

    x_plus = vector<double>(5);
    x_plus(0) = r0;
    x_plus(1) = v0;
    x_plus(2) = 0;
    x_plus(3) = a0;
    x_plus(4) = snr0;

    x_minus = x_plus;

    P_plus = identity_matrix<double>(5) * 5;
    P_minus = identity_matrix<double>(5) * 10;

    t = t0;
    T = 0;

    hasBeenUpdated = false;
    updatesWithoutDetections = 0;

    closest_detection = Detection();

    first_update = true;
}

void Tracker::motionUpdate(double tn, bool updateWihtoutDetection) {
    T = tn-t;
    t_old = tn;

    x_minus = prod(F(), x_plus);

    P_minus = prod(F(), P_plus);
    P_minus = prod(P_minus, trans(F())) + Q();

    if(updateWihtoutDetection) {
        updatesWithoutDetections++;
        hasBeenUpdated = true;
    }
}

void Tracker::informationUpdate(double rn, double an, double vn, double snrn, double dr, double dv) {
    matrix<double> tmp = prod(P_minus, trans(H()));
    matrix<double> S = prod(H(), tmp) + R(snrn, dr, dv);
    matrix<double> K = prod(tmp,matrixInverse(S));

    vector<double> y(4), z(4);
    z(0) = rn; z(1) = vn; z(2) = an; z(3) = snrn;
    y = z - prod(H(),x_minus);
    x_plus = x_minus + prod(K,y);

    matrix<double> tmp2 = prod(K, H());
    P_plus = P_minus - prod(tmp2, P_minus);

    updatesWithoutDetections=0;
    hasBeenUpdated=false;
    t = t_old;

    first_update=false;
}

matrix<double> Tracker::covarianceStateError(double &snrn, double &dr, double &dv, double &dt) {
    // "Motion update"
    matrix<double> tmp = prod(F(dt),P_plus);
    matrix<double> P_m = prod(tmp,trans(F(dt))) + Q(dt);

    // "Information update"
    matrix<double> tmp2 = prod(P_m, trans(H()));
    matrix<double> St = prod(H(), tmp2) + R(snrn, dr, dv);
    matrix<double> Kt = prod(tmp2, matrixInverse(St));
    matrix<double> tmp3 = prod(Kt,H());
    matrix<double> Pt = P_m - prod(tmp3,P_m);

    return Pt;
}

matrix<double> Tracker::covarianceMeasureError(double &snrn, double &dr, double &dv) {
    matrix<double> tmp = prod(H(),P_minus);
    matrix<double> S = prod(tmp,trans(H())) + R(snrn, dr, dv);

    return S;
}


void Tracker::trackPrediction(double &rn, double &an, double &vn, double &snrn) {
    vector<double> v1 = prod(H(), x_minus);
    rn = v1(0);
    vn = v1(1);
    an = v1(2);
    snrn = v1(3);
}

void Tracker::trackEstimate(double &rn, double &an, double &vn, double &snrn) {
    vector<double> v1 = prod(H(), x_plus);
    rn = v1(0);
    vn = v1(1);
    an = v1(2);
    snrn = v1(3);
}

vector<double> Tracker::covariance() {
    vector<double> v1(5);
    for(int i=0; i<5; i++)
        v1(i) = P_minus(i,i);

    return v1;
}

vector<double> Tracker::covariancePlus() {
    vector<double> v1(5);
    for(int i=0; i<5; i++)
        v1(i) = P_plus(i,i);

    return v1;
}

matrix<double> Tracker::F() {
    return F(T);
}

matrix<double> Tracker::F(double &dt) {
    matrix<double> ret = identity_matrix<double>(5);

    ret(0,1) = dt;
    ret(0,2) = (1/std::pow(alpha,2))*(-1+alpha*dt+ std::exp(-alpha*dt));
    ret(1,2) = (1/alpha)*(1-std::exp(-alpha*dt));
    ret(2,2) = std::exp(-alpha*dt);

    return ret;
}

matrix<double> Tracker::H() {
    matrix<double> ret = identity_matrix<double>(4,5);
    ret(2,2) = 0;
    ret(3,4) = 1;

    return ret;
}

matrix<double> Tracker::Q() {
    return Q(T);
}

matrix<double> Tracker::Q(double &dt) {
    matrix<double> q = identity_matrix<double>(5);

//    double at = alpha*dt;
//    q(0,0) = 2*alpha*std::pow(sigmaM,2)*(1/(2*std::pow(alpha,5))) * (1 - std::exp(-2*at) + 2*at + 2*std::pow(at,3)/3 - 2*std::pow(at,2) - 4*at*std::exp(-at));
//    q(0,1) = 2*alpha*std::pow(sigmaM,2)*(1/(2*std::pow(alpha,4))) * (std::exp(-2*at) + 1 - 2*std::exp(-at) + 2*at*std::exp(-at) - 2*at + std::pow(at,2));
//    q(0,2) = 2*alpha*std::pow(sigmaM,2)*(1/(2*std::pow(alpha,3))) * (1 - std::exp(-2*at) - 2*at*std::exp(-at));
//    q(1,0) = 2*alpha*std::pow(sigmaM,2)*(1/(2*std::pow(alpha,4))) * (std::exp(-2*at) + 1 - 2*std::exp(-at) + 2*at*std::exp(-at) - 2*at + std::pow(at,2));
//    q(1,1) = 2*alpha*std::pow(sigmaM,2)*(1/(2*std::pow(alpha,3))) * (4*std::exp(-at) - 3 - std::exp(-2*at) + 2*at);
//    q(1,2) = 2*alpha*std::pow(sigmaM,2)*(1/(2*std::pow(alpha,2))) * (std::exp(-2*at) + 1 - 2*std::exp(-at));
//    q(2,0) = 2*alpha*std::pow(sigmaM,2)*(1/(2*std::pow(alpha,3))) * (1 - std::exp(-2*at) - 2*at*std::exp(-at));
//    q(2,1) = 2*alpha*std::pow(sigmaM,2)*(1/(2*std::pow(alpha,2))) * (std::exp(-2*at) + 1 - 2*std::exp(-at));
//    q(2,2) = 2*alpha*std::pow(sigmaM,2)*(1/(2*alpha))   * (1 - std::exp(-2*at));
    q(0,0) = 2*alpha*std::pow(sigmaM,2)*std::pow(dt,5)/20;
    q(0,1) = 2*alpha*std::pow(sigmaM,2)*std::pow(dt,4)/8;
    q(0,2) = 2*alpha*std::pow(sigmaM,2)*std::pow(dt,3)/6;
    q(1,0) = 2*alpha*std::pow(sigmaM,2)*std::pow(dt,4)/8;
    q(1,1) = 2*alpha*std::pow(sigmaM,2)*std::pow(dt,3)/3;
    q(1,2) = 2*alpha*std::pow(sigmaM,2)*std::pow(dt,2)/2;
    q(2,0) = 2*alpha*std::pow(sigmaM,2)*std::pow(dt,3)/6;
    q(2,1) = 2*alpha*std::pow(sigmaM,2)*std::pow(dt,2)/2;
    q(2,2) = 2*alpha*std::pow(sigmaM,2)*dt;

    q(3,3) = dt;
    q(4,4) = dt;

    return q;
}


matrix<double> Tracker::R(double SNR, double dr, double dv) {
    matrix<double> m2 = identity_matrix<double>(4);
    m2(0,0) = (10*3/2)*std::pow(dr,2)/std::pow(10,SNR/10);
    m2(1,1) = (10*3/std::pow(M_1_PI,2))*std::pow(dv,2)/std::pow(10,SNR/10);
    m2(2,2) = 1/std::pow(10,SNR/10);
    m2(3,3) = 1;

    return m2;
}

matrix<double> Tracker::matrixInverse(matrix<double> &m) {
    if(m.size1()!=m.size2()) {
        std::cerr << "Matrix must be square for inverse" << std::endl;
        return zero_matrix<double>(0);
    }

    size_t matSize = m.size1();

    matrix<double> mInv = identity_matrix<double>(matSize);
    permutation_matrix<size_t> pm(matSize);
    lu_factorize(m, pm);
    lu_substitute(m, pm, mInv);

    return mInv;
}

vector<double> Tracker::distanceInSigmas(double r, double v, double SNR, double dr, double dv) {
    vector<double> ret = zero_vector<double>(2);

    vector<double> z = prod(H(), x_minus);

    double deltaR = std::sqrt(std::pow(z(0) - r, 2));
    double deltaV = std::sqrt(std::pow(z(1) - v, 2));

    matrix<double> Sp = covarianceMeasureError(SNR,dr,dv);

    ret(0) = deltaR / sqrt(Sp(0,0));
    ret(1) = deltaV / sqrt(Sp(1,1));

    return ret;
}

double Tracker::trace(matrix<double>& A)
{
    size_t n = std::min(A.size1(), A.size2());
    double res = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        res += A(i,i);
    }

    return res;
}
