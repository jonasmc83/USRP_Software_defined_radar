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
#ifndef PROCESSING_GPU_H
#define PROCESSING_GPU_H

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include <cufftXt.h>

#define cuPrecisionComplex cuDoubleComplex

void matchedFilterProcessingCUDA_gpu(cuDoubleComplex *signal, cuDoubleComplex *waveform, cuDoubleComplex *window, long int M, long int N);
void dopplerProcessingCUDA_gpu(cuDoubleComplex *signal, long int M, long int N);
void rangeDopplerProcessingCUDA_gpu(cuDoubleComplex *signal, cuDoubleComplex *waveform, cuDoubleComplex *rangeWindow, cuDoubleComplex *dopplerWindow, long int M, long int N);
void rangeDopplerCFARProcessingCUDA_gpu(cuDoubleComplex *signal, cuDoubleComplex *waveform, cuDoubleComplex *rangeWindow, cuDoubleComplex *dopplerWindow, long int M, long int N, int windowLength, int guardInterval);

#endif // PROCESSING_GPU_H
