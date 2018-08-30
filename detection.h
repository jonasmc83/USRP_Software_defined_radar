/*
    USRP_Radar is a software for real time sampling, processing, display and storing
    Copyright (C) 2018  Jonas Myhre Christiansen <jonas-myhre.christiansen@ffi.no>
	
    This file is part of USRP_Radar.

    USRP_Radar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    USRP_Radar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with USRP_Radar.  If not, see <https://www.gnu.org/licenses/>.
*/
	
#ifndef DETECTION_H
#define DETECTION_H

#include <vector>

class Detection {
public:
    Detection();
    void clear();
    void remove_detections(std::vector<long int> sortedIndicesToDelete);

    std::vector<double> range;
    std::vector<double> velocity;
    std::vector<double> snr;
    std::vector<double> azimuth;
    std::vector<double> diff_phase;
    std::vector<long int> rindices;
    std::vector<long int> vindices;
};

#endif // DETECTION_H
