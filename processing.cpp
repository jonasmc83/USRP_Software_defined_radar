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
#include "processing.h"


inline PRECISION windowFunction(long int k, long int N, int windowType=HAMMING_WINDOW) {
    PRECISION w = 0;
    switch(windowType) {
    case HAMMING_WINDOW:
        w = (PRECISION)(0.54 - 0.46*std::cos((2*M_PI*k)/((N)-1)));
        break;
    case BLACKMAN_WINDOW:
        w = (PRECISION)(0.42659 - 0.49656*std::cos((2*M_PI*k)/((N)-1)) + 0.076849*std::cos((4*M_PI*k)/((N)-1)));
        break;
    case BLACKMAN_HARRIS_WINDOW:
        w = (PRECISION)(0.35875 - 0.48829*std::cos((2*M_PI*k)/((N)-1)) + 0.14128*std::cos((4*M_PI*k)/((N)-1)) + 0.01168*std::cos((6*M_PI*k)/((N)-1)));
        break;
    case RECT_WINDOW:
        w = (PRECISION)1;
        break;
    }
    return w;
}

void waveformFFT(std::vector<std::complex<PRECISION>> &Waveform, int N2) {
    kissfft<PRECISION> fft(N2, false);

    std::vector<std::complex<PRECISION>> out(N2);

    fft.transform(&Waveform[0], &out[0]);
    memcpy(&Waveform[0], &out[0], N2 * sizeof(std::complex<PRECISION>));
}

void matchedFilterProcessingKISS(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::complex<PRECISION>> &Waveform, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, size_t N2, size_t Nfull, size_t Noffset, int windowType, int savePulseCompressionData) {
#pragma omp parallel
    {
        kissfft<PRECISION> fftForward((long int)N2, false);
        kissfft<PRECISION> fftReverse((long int)N2, true);

        std::vector<std::complex<PRECISION>> in(N2);
        std::vector<std::complex<PRECISION>> out(N2);

#pragma omp for
        for(long int i=0; i<M; i++) {
            memcpy(&in[0], dataMatrix[i], sizeof(std::complex<PRECISION>)*std::min(N2,Nfull-Noffset));

            if(N2>(Nfull-Noffset)) {
                memset((void*)(&in[0]+Nfull-Noffset), 0, (N2-Nfull+Noffset)*sizeof(std::complex<PRECISION>));
            }

            fftForward.transform(&in[0], &out[0]);

            for(long int k=0; k<N2; k++) {
                out[k] = out[k] * std::conj(Waveform[k]) * windowFunction(k-(int)(N2/2), (int)N2, windowType);
            }

            fftReverse.transform(&out[0], &in[0]);

            memcpy(dataMatrix[i], &in[0], sizeof(std::complex<PRECISION>)*N);

            // If save pulse compressed data
            if(savePulseCompressionData) {
                memcpy(&(saveData[i]->front()), dataMatrix[i], sizeof(std::complex<PRECISION>)*N);
            }
        }
    }
}

void matchedFilterProcessingCUDA(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::complex<PRECISION>> &Waveform, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, size_t N2, size_t Nfull, size_t Noffset, int windowType, int savePulseCompressionData) {

    cuPrecisionComplex *signal = (cuPrecisionComplex*) malloc(sizeof(cuPrecisionComplex)*M*N2);
    cuPrecisionComplex *waveform = (cuPrecisionComplex*) malloc(sizeof(cuPrecisionComplex)*N2);
    cuPrecisionComplex *window = (cuPrecisionComplex*) malloc(sizeof(cuPrecisionComplex)*N2);

    memset(signal, 0, sizeof(cuPrecisionComplex)*M*N2);
    memset(waveform, 0, sizeof(cuPrecisionComplex)*N2);
    for(long int i=0; i<M; i++) {
        memcpy((void*)(&signal[i*N2]), (void*)dataMatrix[i], sizeof(cuPrecisionComplex)*std::min(N2,Nfull-Noffset));
    }
    memcpy((void*) waveform, (void*) (&Waveform[0]), sizeof(cuPrecisionComplex)*N2);
    for(long int i=0; i<N2; i++) {
        window[i].x = windowFunction(i-(int)(N2/2), (int)N2, windowType);
        window[i].y = 0;
    }

    matchedFilterProcessingCUDA_gpu(signal, waveform, window, M, N2);

    for(long int i=0; i<M; i++) {
        memcpy((void*)dataMatrix[i], (void*)(&signal[i*N2]), sizeof(cuPrecisionComplex)*std::min(N2,Nfull-Noffset));

        // If save pulse compressed data
        if(savePulseCompressionData) {
            memcpy(&(saveData[i]->front()), dataMatrix[i], sizeof(std::complex<PRECISION>)*N);
        }
    }

    free(signal);
    free(waveform);
    free(window);
}

void dopplerProcessingKISS(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, int windowType, int saveDopplerData) {
    // Doppler processing
#pragma omp parallel
    {
        kissfft<PRECISION> fftForward((long int)M, false);

        std::vector<std::complex<PRECISION>> in(M);
        std::vector<std::complex<PRECISION>> out(M);

#pragma omp for
        for(long int i=0; i<N; i++) {
            for(long int k=0; k<M; k++) {
                in[k] = dataMatrix[k][i] * windowFunction(k, (long int)M, windowType);
            }

            fftForward.transform(&in[0], &out[0]);

            for(long int k=0; k<M; k++) {
                dataMatrix[k][i] = out[k];
            }
        }

        if(saveDopplerData) {
            for(long int i=0; i<M; i++) {
                memcpy(&(saveData[i]->front()), dataMatrix[i], sizeof(std::complex<PRECISION>)*N);
            }
        }
    }
}

void dopplerProcessingCUDA(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, int windowType, int saveDopplerData) {
    cuPrecisionComplex *in = (cuPrecisionComplex*) malloc(sizeof(cuPrecisionComplex)*M*N);
    for(long int i=0; i<N; i++) {
        for(long int k=0; k<M; k++) {
            in[k+i*M].x = dataMatrix[k][i].real() * windowFunction(k, (long int)M, windowType);
            in[k+i*M].y = dataMatrix[k][i].imag() * windowFunction(k, (long int)M, windowType);
        }
    }
    dopplerProcessingCUDA_gpu(in, M, N);
    for(long int i=0; i<N; i++) {
        for(long int k=0; k<M; k++) {
            dataMatrix[k][i].real(in[k+i*M].x);
            dataMatrix[k][i].imag(in[k+i*M].y);
        }
    }
    free(in);

    if(saveDopplerData) {
        for(long int i=0; i<M; i++) {
            memcpy(&(saveData[i]->front()), dataMatrix[i], sizeof(std::complex<PRECISION>)*N);
        }
    }
}

void FFTProcessingKISS(std::vector<std::complex<PRECISION>> &dataMatrix, size_t M, int windowType) {
    // Doppler processing
//#pragma omp parallel
    {
        kissfft<PRECISION> fftForward((long int)M, false);

        std::vector<std::complex<PRECISION>> in(M);
        std::vector<std::complex<PRECISION>> out(M);

//#pragma omp for
        for(long int k=0; k<M; k++) {
            in[k] = dataMatrix[k] * windowFunction(k, (long int)M, windowType);
        }

        fftForward.transform(&in[0], &out[0]);

        for(size_t k=0; k<M; k++) {
            dataMatrix[k] = out[k];
        }
    }
}

void inverseFFTProcessingKISS(std::vector<std::complex<PRECISION>> &dataMatrix, size_t M, int windowType) {
    // Doppler processing
//#pragma omp parallel
    {
        kissfft<PRECISION> fftReverse((long int)M, true);

        std::vector<std::complex<PRECISION>> in(M);
        std::vector<std::complex<PRECISION>> out(M);

//#pragma omp for
        for(long int k=0; k<M; k++) {
            in[k] = dataMatrix[k] * windowFunction(k, (long int)M, windowType);
        }

        fftReverse.transform(&in[0], &out[0]);

        for(size_t k=0; k<M; k++) {
            dataMatrix[k] = out[k];
        }
    }
}

void rangeDopplerProcessingCUDA(std::vector<std::complex<PRECISION>*> &dataMatrix, std::vector<std::complex<PRECISION>> &Waveform, std::vector<std::vector<std::complex<PRECISION>>*> &saveData, size_t M, size_t N, size_t N2, size_t Nfull, size_t Noffset, int rangeWindowType, int dopplerWindowType, int saveDataFlag, int CFARFlag, int windowLength, int guardInterval) {
    cuPrecisionComplex *signal = (cuPrecisionComplex*) malloc(sizeof(cuPrecisionComplex)*M*N2);
    cuPrecisionComplex *waveform = (cuPrecisionComplex*) malloc(sizeof(cuPrecisionComplex)*N2);
    cuPrecisionComplex *rangeWindow = (cuPrecisionComplex*) malloc(sizeof(cuPrecisionComplex)*N2);
    cuPrecisionComplex *dopplerWindow = (cuPrecisionComplex*) malloc(sizeof(cuPrecisionComplex)*M);

    memset(signal, 0, sizeof(cuPrecisionComplex)*M*N2);
    memset(waveform, 0, sizeof(cuPrecisionComplex)*N2);
    for(long int i=0; i<M; i++) {
        memcpy((void*)(&signal[i*N2]), (void*)dataMatrix[i], sizeof(cuPrecisionComplex)*std::min(N2,Nfull-Noffset));

        dopplerWindow[i].x = windowFunction(i, (long int)M, dopplerWindowType);
        dopplerWindow[i].y = 0;
    }
    memcpy((void*) waveform, (void*) (&Waveform[0]), sizeof(cuPrecisionComplex)*N2);
    for(long int i=0; i<N2; i++) {
        rangeWindow[i].x = windowFunction(i-(int)(N2/2), (int)N2, rangeWindowType);
        rangeWindow[i].y = 0;
    }

    if(CFARFlag) {
        rangeDopplerCFARProcessingCUDA_gpu(signal, waveform, rangeWindow, dopplerWindow, M, N2, windowLength, guardInterval);
    } else {
        rangeDopplerProcessingCUDA_gpu(signal, waveform, rangeWindow, dopplerWindow, M, N2);
    }

    for(long int i=0; i<M; i++) {
        memcpy((void*)dataMatrix[i], (void*)(&signal[i*N2]), sizeof(cuPrecisionComplex)*std::min(N2,Nfull-Noffset));

        // If save data
        if(saveDataFlag) {
            memcpy(&(saveData[i]->front()), dataMatrix[i], sizeof(std::complex<PRECISION>)*N);
        }
    }

    free(signal);
    free(waveform);
    free(rangeWindow);
    free(dopplerWindow);
}

void cfarProcessing(std::vector<std::complex<PRECISION>*> &dataMatrix, long int M, long int N, long int windowLength, long int guardInterval)  {
#pragma omp parallel

#pragma omp for
    for(size_t k=0; k<M; k++) {
        double Sum=0;
        std::vector<PRECISION> *Cfar = new std::vector<PRECISION>(N,std::nan(""));
        for(size_t i=guardInterval; i<N; i++) {
            Sum += std::pow(std::abs(dataMatrix[k][i-guardInterval]),2);
            if(i>guardInterval+windowLength) {
                Sum -= std::pow(std::abs(dataMatrix[k][i-(guardInterval+windowLength+1)]),2);
            }
            (*Cfar)[i] = Sum;
        }
        for(size_t i=1; i<N; i++) {
            if(i<guardInterval+windowLength) {
                dataMatrix[k][i] /= std::sqrt(2*(*Cfar)[i]/i);
            } else {
                dataMatrix[k][i] /= std::sqrt(2*(*Cfar)[i]/windowLength);
            }
        }
        if(Cfar!=NULL) {
            delete Cfar;
            Cfar=NULL;
        }
    }
}

int upsampleFT(std::vector<PRECISION> &profile_out, std::vector<PRECISION> &profile_in, int UpsamplingFactor) {
    size_t Nin = profile_in.size();
    size_t Nout = Nin*UpsamplingFactor;

    size_t Nin2f = std::floor(Nin/2);
    size_t Nin2c = std::ceil(Nin/2);

    if(profile_in.size()<2) {
        std::cerr << "Upsampling vector must contain more than one element" << std::endl;
        return -1;
    }
    if(Nout != profile_out.size()) {
        std::cerr << "profile_out must contain Nin*UpsamplingFactor elements" << std::endl;
        return -1;
    }

    std::vector<std::complex<PRECISION>> tmp_in(Nin), tmp_in2(Nin);
    std::vector<std::complex<PRECISION>> tmp_out2(Nout), tmp_out(Nout, std::complex<PRECISION>(0,0));

    for(long int k=0; k<Nin; k++) {
        tmp_in[k] = std::complex<PRECISION>(profile_in[k],0);
    }

    kissfft<PRECISION> fftForward((long int)Nin, false);
    fftForward.transform(&tmp_in[0], &tmp_in2[0]);

    memcpy(&tmp_out[0], &tmp_in2[0], sizeof(std::complex<PRECISION>)*Nin2f);
    memcpy(&tmp_out[Nout-Nin2c], &tmp_in2[Nin2f+1], sizeof(std::complex<PRECISION>)*Nin2c);

    kissfft<PRECISION> fftReverse((long int)Nout, true);
    fftReverse.transform(&tmp_out[0], &tmp_out2[0]);

    for(long int k=0; k<Nout; k++) {
        profile_out[k] = std::abs(tmp_out2[k]);
    }

    return 1;
}

double absval(double x) {
    return std::sqrt(std::pow(x,2));
}
