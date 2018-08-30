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
#ifndef WAVEFORM_H
#define WAVEFORM_H

#include <complex>
#include <vector>
#include <fstream>

#include "precision.h"

#define WAVEFORM_ERROR_VECTOR_TOO_SHORT -1

#define WAVEFORM_TYPE_LFM 0
#define WAVEFORM_TYPE_RECT 1
#define WAVEFORM_TYPE_FROM_FILE 2

class Waveform
{
public:
    Waveform();

    double Bandwidth;
    double f0;
    double fs;
    double tau;
    int Waveform_type;
    int Waveform_mask;
    std::string fileName;

    long int generateWaveform(std::vector<std::complex<PRECISION>> *wfVector);
};

#endif // WAVEFORM_H
