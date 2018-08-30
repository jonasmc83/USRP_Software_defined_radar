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
#include "waveform.h"

Waveform::Waveform()
{
    Waveform_type = WAVEFORM_TYPE_LFM;
}

long int Waveform::generateWaveform(std::vector<std::complex<PRECISION>> *wfVector) {
    long int Nwaveform = (long int) (tau * fs);
    long int k;

    std::vector<double> phi(Nwaveform,0);
    std::vector<PRECISION> *from_file = NULL;

    if (wfVector->size()<Nwaveform) {
        return WAVEFORM_ERROR_VECTOR_TOO_SHORT;
    }

    if(Waveform_type==WAVEFORM_TYPE_FROM_FILE) {
        std::ifstream in(fileName, std::ios::binary);

        in.seekg (0, in.end);
        long int length = in.tellg();
        in.seekg (0, in.beg);

        from_file = new std::vector<PRECISION>(length/4);
        in.read((char*) &(from_file->front()), length);
        in.close();

        Nwaveform = (size_t) length/8;
    }

    switch(Waveform_type) {
    // Quadratic phase evolution (linear freq)
    case WAVEFORM_TYPE_LFM:
        for(k=1;k<Nwaveform;k++) {
            phi[k] = phi[k-1] + (2*M_PI)*((k-Nwaveform/2)*Bandwidth/Nwaveform)/fs;
        }
        break;
    // Const phase 0, zero frequency pules
    case WAVEFORM_TYPE_RECT:
        break;
    }

    for(k=0; k<Nwaveform; k++) {
        // If from file, storing binary waveform in wfVector
        if(Waveform_type==WAVEFORM_TYPE_FROM_FILE) {
            (*wfVector)[k] = std::complex<PRECISION>((*from_file)[2*k], (*from_file)[2*k+1]);
        } else {
            (*wfVector)[k] = std::exp(std::complex<PRECISION>(0,phi[k]));
        }
    }

    if(Waveform_type==WAVEFORM_TYPE_FROM_FILE && from_file!=NULL) {
        delete from_file;
    }

    return Nwaveform;
}

