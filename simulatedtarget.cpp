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
#include "simulatedtarget.h"

#include <algorithm>
#include <complex>

SimulatedTarget::SimulatedTarget(double x0, double y0, double vx, double vy, double RCS)
{
    this->x0 = x0;
    this->y0 = y0;
    this->vx = vx;
    this->vy = vy;
    this->x = x0;
    this->y = y0;
    this->rcs = RCS;
}

void SimulatedTarget::updateTarget(double T) {
    x += vx*T;
    y += vy*T;
}

double SimulatedTarget::getRange() {
    return std::sqrt( std::pow(x,2) + std::pow(y,2) );
}

double SimulatedTarget::getVelocity() {
    double delta = M_PI/4 - std::atan2(vy,vx);
    return std::sqrt( std::pow(vx,2) + std::pow(vy,2) )*std::cos(delta - getAzimuth());
}

double SimulatedTarget::getAzimuth() {
    return M_PI/4 - std::atan2(y,x);
}

double SimulatedTarget::getRCS() {
    return rcs;
}
