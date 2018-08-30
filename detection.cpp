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
	
#include "detection.h"

Detection::Detection()
{

}


void Detection::clear() {
    range.clear();
    velocity.clear();
    snr.clear();
    azimuth.clear();
    diff_phase.clear();
    rindices.clear();
    vindices.clear();
}

void Detection::remove_detections(std::vector<long int> sortedIndicesToDelete) {
    std::vector<double> tmpRange;
    std::vector<double> tmpVelocity;
    std::vector<double> tmpSNR;
    std::vector<double> tmpAzimuth;
    std::vector<double> tmpDiff;
    std::vector<long int> tmpRindices;
    std::vector<long int> tmpVindices;

    for(long int i=range.size()-1; i>=0; i--) {
        if(i!=sortedIndicesToDelete.back()) {
            tmpRange.push_back(range[i]);
            tmpVelocity.push_back(velocity[i]);
            tmpSNR.push_back(snr[i]);
            tmpAzimuth.push_back(azimuth[i]);
            tmpDiff.push_back(diff_phase[i]);
            tmpRindices.push_back(rindices[i]);
            tmpVindices.push_back(vindices[i]);
        } else {
            sortedIndicesToDelete.pop_back();
        }
    }

    range.clear(); range = tmpRange;
    velocity.clear(); velocity = tmpVelocity;
    snr.clear(); snr = tmpSNR;
    azimuth.clear(); azimuth = tmpAzimuth;
    diff_phase.clear(); diff_phase = tmpDiff;
    rindices.clear(); rindices = tmpRindices;
    vindices.clear(); vindices = tmpVindices;
}
