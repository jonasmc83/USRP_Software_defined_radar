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
#ifndef SIMULATEDTARGET_H
#define SIMULATEDTARGET_H

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

class SimulatedTarget
{
private:
    double x0,y0,x,y,vx,vy,rcs;

public:
    SimulatedTarget(double x0, double y0, double vx, double vy, double rcs);
    void updateTarget(double T);
    double getRange();
    double getVelocity();
    double getAzimuth();
    double getRCS();
};

#endif // SIMULATEDTARGET_H
