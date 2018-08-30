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
#ifndef PROCESSING_H
#define PROCESSING_H

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define HAMMING_WINDOW 0
#define BLACKMAN_WINDOW 1
#define BLACKMAN_HARRIS_WINDOW 2
#define RECT_WINDOW 3

#include <omp.h>

#include <iostream>
#include <vector>
#include <cstring>
#include <complex>
#include <algorithm>
#include <kissfft.hh>
#include <cmath>

#include "precision.h"

#ifdef FLOAT
    #include "processing_gpu.h"
#else
    #include "processing_gpu_double.h"
#endif

inline PRECISION windowFunction(long int k, long int N, int windowType);
void waveformFFT(std::vector<std::complex<PRECISION>> &Waveform, int N2);

void matchedFilterProcessingKISS(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::complex<PRECISION>> &Waveform, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, size_t N2, size_t Nfull, size_t Noffset, int windowType, int savePulseCompressionData);
void matchedFilterProcessingCUDA(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::complex<PRECISION>> &Waveform, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, size_t N2, size_t Nfull, size_t Noffset, int windowType, int savePulseCompressionData);

void dopplerProcessingKISS(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, int windowType, int saveDopplerData);
void dopplerProcessingCUDA(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, int windowType, int saveDopplerData);

void FFTProcessingKISS(std::vector<std::complex<PRECISION>> &dataMatrix, size_t M, int windowType);
void inverseFFTProcessingKISS(std::vector<std::complex<PRECISION>> &dataMatrix, size_t M, int windowType);

void rangeDopplerProcessingCUDA(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::complex<PRECISION>> &Waveform, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, size_t N2, size_t Nfull, size_t Noffset, int rangeWindowType, int dopplerWindowType, int saveDataFlag, int CFARFlag, int windowLength=0, int guardInterval=0);

void cfarProcessing(std::vector<std::complex<PRECISION>*> &dataMatrix, long int M, long int N, long int windowLength, long int guardInterval);

int upsampleFT(std::vector<PRECISION> &profile_out, std::vector<PRECISION> &profile_in, int UpsamplingFactor);

double absval(double x);

#endif // PROCESSING_H
