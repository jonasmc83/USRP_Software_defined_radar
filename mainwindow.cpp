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
#include "mainwindow.h"
#include "ui_mainwindow.h"

/***********************************************************************
 * transmit_pulse_train_worker function
 * A function to be used as a boost::thread_group thread for transmitting
 **********************************************************************/
void transmit_pulse_train_worker(uhd::usrp::multi_usrp::sptr usrp, std::vector<std::complex<PRECISION>*> buff_ptrs, size_t N, size_t M, uhd::time_spec_t time_to_start){

    //create a transmit streamer
    static bool first = true;

    static uhd::stream_args_t tx_stream_args(USRP_PRECISION, "sc16");
    if(first)
        tx_stream_args.channels.push_back(0);

    static uhd::tx_streamer::sptr tx_stream;
    if(first)
        tx_stream = usrp->get_tx_stream(tx_stream_args);

    // create metadata structure
    static uhd::tx_metadata_t tx_md;
    tx_md.start_of_burst = true;
    tx_md.end_of_burst = false;
    tx_md.has_time_spec = true;
    tx_md.time_spec = time_to_start;

    // transmit loop (send pulses back to back)
    for(long int k = 0; k < M; k++) {
        //send the entire contents of the buffer
        tx_stream->send(buff_ptrs, N, tx_md, 0.5);

        tx_md.start_of_burst = false;
        tx_md.has_time_spec = false;
    }

    //send end of burst packet
    tx_md.has_time_spec = false;
    tx_md.end_of_burst = true;
    tx_stream->send("", 0, tx_md);

    first=false;
}

/***********************************************************************
 * receive_pulse_train function
 **********************************************************************/
long int receive_pulse_train(uhd::usrp::multi_usrp::sptr usrp, std::vector<std::complex<PRECISION>*> buff_ptrs, size_t len, uhd::time_spec_t time_to_start){
    // Setup stream args
    static bool first = true;

    static size_t channels = usrp->get_rx_num_channels();
    static std::vector<size_t> channel_vector;
    static uhd::stream_args_t stream_args(USRP_PRECISION,"sc16");

    if(first) {
        if(channel_vector.size()==0){
            for(size_t i=0; i<channels; i++)
                channel_vector.push_back(i);
        }
        stream_args.channels = channel_vector;
    }

    static uhd::rx_streamer::sptr rx_stream;

    if(first)
        rx_stream= usrp->get_rx_stream(stream_args);

    // creating rx metadata structure
    uhd::rx_metadata_t md;

    //setup streaming
    uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
    stream_cmd.num_samps = len;
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = time_to_start;
    rx_stream->issue_stream_cmd(stream_cmd);

    // Sampling
    size_t max_samps = rx_stream->get_max_num_samps();

    size_t num_acc_samps = 0; //number of accumulated samples
    std::vector<std::complex<PRECISION> *> buff_ptrs2(buff_ptrs.size());
    while(num_acc_samps < len){
        // Moving storing pointer to correct location
        for (size_t i = 0; i < channels; i++) buff_ptrs2[i] = &(buff_ptrs[i][num_acc_samps]);

        // sampling data
        size_t samps_to_recv = std::min(len - num_acc_samps, max_samps);
        size_t num_rx_samps = rx_stream->recv(buff_ptrs2, samps_to_recv, md, 0.5);

        num_acc_samps += num_rx_samps;

        //handle the error code
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) break;
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE) break;
    }

    first = false;

    if(num_acc_samps<len)
        return -((long int)num_acc_samps);
    else
        return (long int)num_acc_samps;
}

// Main processing thread running
void pulse_doppler_worker(MainWindow *win) {
    boost::thread_group transmit_thread;

    size_t CPINum = 0;

    double pulse_length = win->ui->pulseLength->text().toDouble();
    double tx_rate = win->ui->txRate->text().toDouble();
    double rx_rate = win->ui->rxRate->text().toDouble();
    double PRI = win->ui->PRI->text().toDouble();
    double f0 = win->ui->carrierFrequency->text().toDouble();
    double dr = (3e8/(2*rx_rate));
    double rxTimeGPS;
    size_t M = win->ui->numPulses->text().toULong();
    size_t Nfull = (size_t) (PRI*rx_rate);
    size_t Noffset = (size_t) (win->ui->rangeOffset->text().toDouble()/dr);
    size_t Moffset = (size_t) (win->ui->timeOffset->text().toDouble()/PRI);
    size_t N = std::min((size_t)(win->ui->maxRange->text().toDouble()/(3e8/(2*rx_rate))), Nfull) - Noffset;
    size_t Nwaveform = (size_t) (pulse_length * tx_rate);
    size_t len = (size_t) (M*Nfull);
    size_t N2 = std::pow(2, std::ceil(std::log2(std::max(N,Nwaveform)+Nwaveform)));
    double BW = win->ui->BW->text().toDouble() / tx_rate;
    double dv = 3e8/(2*PRI*f0*M);
    bool saveData = false;

    double MinUpdateRate = 0.3;
    double WantedUpdateRate;

    double wantedPRI;

    Waveform wf;

    std::vector<std::complex<PRECISION>> *rx_buff_ch1=NULL, *rx_buff_ch2=NULL, *tx_buff=NULL;

    if(M>MAX_NUM_PULSES) {
        std::cout << "Number of pulses is too large: MAX=" << MAX_NUM_PULSES << std::endl;
        QMetaObject::invokeMethod(win, SLOT(on_pushButton_2_clicked));
        return;
    }

    // Pointer buffers
    std::vector<std::complex<PRECISION>*> rx_buff_ptrs;
    std::vector<std::complex<PRECISION>*> tx_buff_ptrs;

    std::vector<std::complex<PRECISION>*> dataMatrixCh1(MAX_NUM_PULSES);
    std::vector<std::complex<PRECISION>*> dataMatrixCh2(MAX_NUM_PULSES);

    // Processing loop
    uhd::time_spec_t time_to_start;
    boost::chrono::system_clock::time_point last = boost::chrono::system_clock::now();
    while(win->running) {
        long int num_samples_received;

        boost::chrono::system_clock::time_point now = boost::chrono::system_clock::now();
        boost::chrono::duration<double> sec = now - last;
        last = now;

        std::vector<double> updateRateVector;
        std::vector<double> costFunctionVector;

        // Check if we will save data
        if(win->ui->saveDataToFile->isChecked())
            saveData=true;
        else
            saveData=false;

        // Track motion update
        boost::chrono::duration<double> t = boost::chrono::system_clock::now() - win->start;
        for(long int i=0; i<win->trackers.size(); i++) {
            win->trackers[i]->motionUpdate(t.count(), true);
            double r,a,v,s;
            win->trackers[i]->trackPrediction(r,a,v,s);
            std::cout << "Track #" << i << " prediction at " << t.count() << "s : " << r << " m, " << v << " m/s, " << a << " rad, " << s << " dB" << std::endl;
        }

        if(CPINum>0)
            std::cout << "Update rate: " << sec.count() << " seconds at #" << CPINum << " CPI" << std::endl;

        { // Updating waveform parameters
            WantedUpdateRate = MinUpdateRate;
            if(win->ui->checkUpdateInterval->isChecked()) {
                WantedUpdateRate = win->ui->txtUpdateInterval->text().toDouble();
                std::cout << "Wanted update rate set to " << WantedUpdateRate << "s" << std::endl;
            }

            wantedPRI = win->ui->PRI->text().toDouble();
            wantedPRI = std::min(std::max(wantedPRI, (double)(MIN_PRI)), (double)(MAX_PRI));

            pulse_length = win->ui->pulseLength->text().toDouble();
            tx_rate = win->ui->txRate->text().toDouble();
            rx_rate = win->ui->rxRate->text().toDouble();
            PRI = wantedPRI;
            win->ui->PRI->text().toDouble();
            f0 = win->ui->carrierFrequency->text().toDouble();
            M = win->ui->numPulses->text().toULong();
            Nfull = (size_t) (PRI*rx_rate);
            N = std::min((size_t)(win->ui->maxRange->text().toDouble()/(3e8/(2*rx_rate))), Nfull-Noffset);
            Nwaveform = (size_t) (pulse_length * tx_rate);
            len = (size_t) ((M+Moffset)*Nfull);
            N2 = std::pow(2, std::ceil(std::log2(std::max(N,Nwaveform)+Nwaveform)+1));
            BW = win->ui->BW->text().toDouble() / tx_rate; // digital bw

            // waveform object stuff
            wf.Bandwidth = win->ui->BW->text().toDouble();
            wf.f0 = win->ui->carrierFrequency->text().toDouble();
            wf.fs = win->ui->txRate->text().toDouble();
            wf.tau = win->ui->pulseLength->text().toDouble();
            wf.Waveform_type = win->ui->waveformType->currentIndex();

            // Checking if CPI is too long or short
            if(wantedPRI*M<MIN_CPI) {
                M = std::ceil(MIN_CPI/wantedPRI);
                std::cout << "CPI is too low, changing num pulses to: " << M << std::endl;
            } else if(wantedPRI*M>MAX_CPI) {
                M = std::floor(MAX_CPI/wantedPRI);
                std::cout << "CPI is too high, changing num pulses to: " << M << std::endl;
            }

            // Calculating range and velocity resolution
            dr = (3e8/(2*rx_rate));
            dv = 3e8/(2*PRI*f0*M);

            // Storing waveform parameters
            win->prfV.push_back(1/PRI);
            win->bwV.push_back(BW);
            win->plV.push_back(pulse_length);
            win->npV.push_back(M);
            win->tV.push_back(t.count());

            // Displaying unambiguous range and velocity
            double Runamb = (3e8 * PRI)/2;
            double Vunamb = 3e8 / (2*PRI*f0*2);
            win->ui->unambRange->setText(QString::number(Runamb));
            win->ui->unambVelocity->setText(QString::number(Vunamb));
        }

        // Allocating buffers
        if(rx_buff_ch1!=NULL || rx_buff_ch2!=NULL || tx_buff!=NULL) {
            std::cout << "buffers not deleted" << std::endl;
            return;
        }

        // Allocating rx and tx buffers
        rx_buff_ch1 = new std::vector<std::complex<PRECISION>>(len, std::complex<PRECISION>(0,0));
        rx_buff_ch2 = new std::vector<std::complex<PRECISION>>(len, std::complex<PRECISION>(0,0));
        tx_buff = new std::vector<std::complex<PRECISION>>(Nfull, std::complex<PRECISION>(0,0));

        std::vector<std::complex<PRECISION>> vecWaveform(N2, std::complex<PRECISION>(0,0));

        { // Create waveform
            memset(&(tx_buff->front()), 0, sizeof(std::complex<PRECISION>)*Nfull);

            // if waveform from file, store filname to waveform object
            if(wf.Waveform_type==WAVEFORM_TYPE_FROM_FILE) {
                wf.fileName = win->ui->txtFileName->text().toStdString();
            }

            // Generate wavform
            long int retVal = wf.generateWaveform(tx_buff);
            if(retVal<0) {
                std::cerr << "waveform generation error" << std::endl;
                break;
            }
            std::copy(tx_buff->begin(), tx_buff->begin()+retVal, vecWaveform.begin());

            // Waveform FFT
            waveformFFT(vecWaveform, N2);


            if(win->waveform!=NULL) {
                std::cout << "waveform not deleted!" << std::endl;
                return;
            }
            win->waveform = new std::vector<std::complex<PRECISION>>(Nfull);
            win->waveformDim = Nfull;

            // Storing waveform temporal and spectrum for display
            if(win->waveformSpectrum!=NULL) {
                std::cout << "waveformSpectrum no deleted!" << std::endl;
                return;
            }
            win->waveformSpectrum = new std::vector<std::complex<PRECISION>>(N2);
            win->waveformSpectrumDim = N2;
            std::copy(tx_buff->begin(), tx_buff->begin()+Nfull, win->waveform->begin());
            std::copy(vecWaveform.begin(), vecWaveform.begin()+N2, win->waveformSpectrum->begin());
        }

        // Setting up pointer vectors
        rx_buff_ptrs.clear();
        tx_buff_ptrs.clear();
        rx_buff_ptrs.push_back(&rx_buff_ch1->front());
        rx_buff_ptrs.push_back(&rx_buff_ch2->front());
        tx_buff_ptrs.push_back(&tx_buff->front());

        // Setting up data, either simulation or sampling
        std::cout << "Sending-receiving data";
        boost::chrono::system_clock::time_point rxtx_start = boost::chrono::system_clock::now();

        // Single target simulation (could be buggy)
        if(win->simulation) {
            std::mt19937 generator(std::random_device{}());
            auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0},
                                  std::mt19937(std::random_device{}()));

            std::vector<double> Swerling;
            for(long int l=0; l<win->targets.size(); l++) {
                win->targets[l]->updateTarget(sec.count()+0.01*dist());
                // Swerling1 (scan to scan)
                Swerling.push_back(std::abs(dist()));
                std::cout << "Simulated target at range " << win->targets[l]->getRange() << " m and velocity " << win->targets[l]->getVelocity() << " m/s" << std::endl;
            }

            double roffset = dist()*0.1;
            double voffset = dist()*0.1;

            for(size_t i=0; i<(M+Moffset); i++) {
                for(size_t k=0; k<N; k++) {
                    // Adding thermal noise
                    (*rx_buff_ch1)[i*Nfull+k] = std::complex<PRECISION>(5e-5*dist(), 5e-5*dist());
                    (*rx_buff_ch2)[i*Nfull+k] = std::complex<PRECISION>(5e-5*dist(), 5e-5*dist());

                    for(int l=0; l<win->targets.size(); l++) {
                        long int Kstart = (long int) ((win->targets[l]->getRange()+roffset)/dr);
                        if(k>=Kstart && k<Kstart+Nwaveform) {
                            double fd = (win->targets[l]->getVelocity()+voffset)* 2 / (3e8/f0);
                            std::complex<PRECISION> SNR1(Swerling[l]*win->targets[l]->getRCS()*1e3/std::pow(win->targets[l]->getRange(),4),0);
                            std::complex<PRECISION> SNR2 = SNR1;

                            (*rx_buff_ch1)[i*Nfull+k] += (*tx_buff)[k-Kstart] * std::exp(std::complex<PRECISION>(0,((-2*M_PI) * (i*PRI) * fd))) * std::sqrt(SNR1);
                            (*rx_buff_ch2)[i*Nfull+k] += (*tx_buff)[k-Kstart] * std::exp(std::complex<PRECISION>(0,((-2*M_PI) * (i*PRI) * fd))) * std::sqrt(SNR2);
                        }
                    }
                }
            }

            num_samples_received=(long int)len;
        } else { // USRP transmit and receive pulse trains
            // Saving time instant
            double timeSendTX = ((double)M*N*sizeof(PRECISION)*16)/20e9; // Number of bits divided by 20Gbs

            if(CPINum>0) {
                uhd::time_spec_t t_diff = win->usrp->get_time_now() - time_to_start;
                if(WantedUpdateRate>0) {
                    if(t_diff.get_real_secs()>WantedUpdateRate-timeSendTX) {
                        std::cout << "Update rate larger than wanted update rate from tracker, adjust goals" << std::endl;
                    } else {
                        while(t_diff.get_real_secs()<WantedUpdateRate-timeSendTX) {
                            boost::this_thread::sleep(boost::posix_time::milliseconds(5));
                            t_diff = win->usrp->get_time_now() - time_to_start;
                        }
                    }
                }
                vector<double> track_cov = zero_vector<double>(5);
                double r=0,a=0,v=0,s=0;
                if(win->trackers.size()>0) {
                    track_cov = win->trackers[0]->covariancePlus();
                    win->trackers[0]->trackPrediction(r,a,v,s);
                }
            }
            time_to_start = win->usrp->get_time_now() + uhd::time_spec_t(timeSendTX+0.05);
            rxTimeGPS = time_to_start.get_real_secs();

            // Start transmit thread
            transmit_thread.create_thread(boost::bind(&transmit_pulse_train_worker, win->usrp, tx_buff_ptrs, Nfull, M+Moffset, time_to_start));

            // Receive data function
            num_samples_received = receive_pulse_train(win->usrp, rx_buff_ptrs,len,time_to_start);

            // Close transmit thread
            transmit_thread.join_all();
        }

        boost::chrono::system_clock::time_point rxtx_stop = boost::chrono::system_clock::now();
        boost::chrono::duration<double> rxtx_dur = rxtx_stop - rxtx_start;

        std::cout << " - process lasted " << rxtx_dur.count() << " seconds" << std::endl;

        if(num_samples_received<0) {
            std::cout  << "Error occured during sampling" << std::endl;
            std::cout << "Received " << -num_samples_received << " samples out of " << len << std::endl;
        } else {
            if(win->simulation)
                Noffset=0;

            std::vector<std::vector<std::complex<PRECISION>>*> saveDataCh1(M);
            std::vector<std::vector<std::complex<PRECISION>>*> saveDataCh2(M);

            for(size_t i=0; i<M; i++) {
                dataMatrixCh1[i] = &((*rx_buff_ch1)[(i+Moffset)*Nfull+Noffset]);
                dataMatrixCh2[i] = &((*rx_buff_ch2)[(i+Moffset)*Nfull+Noffset]);

                if(saveData) {
                    saveDataCh1[i] = new std::vector<std::complex<PRECISION>>(N,std::complex<PRECISION>(0,0));
                    saveDataCh2[i] = new std::vector<std::complex<PRECISION>>(N,std::complex<PRECISION>(0,0));
                }

                if(win->ui->saveDataType->currentIndex() == SAVE_TYPE_RAW && saveData) {
                    memcpy(&(saveDataCh1[i]->front()), dataMatrixCh1[i], sizeof(std::complex<PRECISION>)*N);
                    memcpy(&(saveDataCh2[i]->front()), dataMatrixCh2[i], sizeof(std::complex<PRECISION>)*N);
                }
            }

            // range-Doppler and CFAR processing
            {
                std::cout << "Starting range-Doppler processing";
                boost::chrono::system_clock::time_point pc_start = boost::chrono::system_clock::now();

                int rangeWindowType = win->ui->windowTypeRange->currentIndex();
                int dopplerWindowType = win->ui->windowTypeDoppler->currentIndex();
                int windowLength = win->ui->windowLength->text().toInt();
                int guardInterval = win->ui->guardInterval->text().toInt();
                int CFARFlag = win->ui->CFAR->isChecked();
                int saveDataFlag = win->ui->saveDataType->currentIndex() == SAVE_TYPE_RANGE_DOPPLER && saveData;

                // Performing processing on GPU (CUDA), change for CPU
                rangeDopplerProcessingCUDA(dataMatrixCh1, vecWaveform, saveDataCh1, M, N, N2, Nfull, Noffset, rangeWindowType, dopplerWindowType, saveDataFlag, CFARFlag, windowLength, guardInterval);
                rangeDopplerProcessingCUDA(dataMatrixCh2, vecWaveform, saveDataCh2, M, N, N2, Nfull, Noffset, rangeWindowType, dopplerWindowType, saveDataFlag, CFARFlag, windowLength, guardInterval);

                boost::chrono::system_clock::time_point pc_stop = boost::chrono::system_clock::now();
                boost::chrono::duration<double> pc_dur = pc_stop - pc_start;

                std::cout << " - processing lasted " << pc_dur.count() << " seconds" << std::endl;
            }

            // Detection
            {
                std::cout << "Starting detection processing";
                boost::chrono::system_clock::time_point det_time_start = boost::chrono::system_clock::now();

                double Threshold = win->ui->threshold->text().toDouble();

                double vWFSize = dv*2;
                double rWFSize = dr*2;

                long int det_start = (long int)(pulse_length*rx_rate);
                if(win->ui->cwCheck->isChecked())
                    det_start=0;

                double NoiseCh1 = 0;
                double NoiseCh2 = 0;

                if(!win->ui->CFAR->isChecked()) {
                    long int Nsum = N-det_start;
                    for(long int k=det_start; k<N; k++) {
                        NoiseCh1 += std::norm(dataMatrixCh1[(long int)(M/2-M/8)][k]);
                        NoiseCh2 += std::norm(dataMatrixCh2[(long int)(M/2-M/8)][k]);
                    }
                    if(NoiseCh1>0)
                        NoiseCh1 = 10*std::log10(2*NoiseCh1/Nsum);
                    else
                        NoiseCh1 = 100;

                    if(NoiseCh2>0)
                        NoiseCh2 = 10*std::log10(2*NoiseCh2/Nsum);
                    else
                        NoiseCh2 = 100;
                }
                win->NoiseFloorCh1 = NoiseCh1;
                win->NoiseFloorCh2 = NoiseCh2;

                win->Detections->clear();

                Detection det;
                int minVelInd = (int) std::min((long int)(win->ui->minVelocity->text().toDouble()/dv), (long int)(M/2));
                for(long int i=minVelInd+1; i<M-minVelInd; i++) { // do not do detection at zero doppler if minVelInd>0

                    for(long int k=det_start; k<(int)N; k++) {
                        double temp, Noise;
                        if(win->ui->channelPlot->currentIndex() == PLOT_CHANNEL_1) {
                            temp = 20*std::log10(std::abs(dataMatrixCh1[i][k]));
                            Noise = NoiseCh1;
                        } else {
                            temp = 20*std::log10(std::abs(dataMatrixCh2[i][k]));
                            Noise = NoiseCh2;
                        }
                        if(temp > (Noise+Threshold)) {
                            double ii = (double)i;
                            if(i>(int)(M/2))
                                ii -= (double)M;
                            double diff_phase = std::arg(dataMatrixCh1[i][k] * std::conj(dataMatrixCh2[i][k]));

                            det.range.push_back((double) (k*dr));
                            det.velocity.push_back((double)(-ii*dv));
                            det.snr.push_back(temp-Noise);
                            det.azimuth.push_back(win->ppiWindow->phaseToAzimuth(diff_phase, f0));
                            det.diff_phase.push_back(diff_phase);
                            det.vindices.push_back(i);
                            det.rindices.push_back(k);
                        }
                    }
                }

                // Allocating variables for range and velocity increased accuracy via FFT interpolation
                long int UpsamplePoints = 7;
                long int UpsampleFactor = 100;
                long int UP2 = std::floor(UpsamplePoints/2);
                std::vector<PRECISION> rin(UpsamplePoints,0), rout(UpsamplePoints*UpsampleFactor,0);
                std::vector<PRECISION> vin(UpsamplePoints,0), vout(UpsamplePoints*UpsampleFactor,0);

                double NCh1Lin = std::pow(10, NoiseCh1/20);
                double NCh2Lin = std::pow(10, NoiseCh2/20);

                if(det.snr.size()>100000) {
                    det.clear();
                    std::cout << "Too many detections! Aborting processing for dwell" << std::endl;
                }

                while(det.snr.size()>0) {
                    std::vector<double>::iterator result = std::max_element(det.snr.begin(), det.snr.end());
                    long int SNRMaxInd = std::distance(det.snr.begin(), result);

                    // Upsampling to find accurate estimate of range and velocity
                    long int m=det.vindices[SNRMaxInd];
                    long int n=det.rindices[SNRMaxInd];

                    for(long int k=-UP2; k<=UP2; k++) {
                        if(win->ui->channelPlot->currentIndex() == PLOT_CHANNEL_1) {
                            if(n+k>=0 && n+k<Nfull) {
                                rin[k+UP2] = std::abs(dataMatrixCh1[m][n+k])/NCh1Lin;
                            } else {
                                rin[k+UP2] = 0;
                            }
                            if(m+k>=0 && m+k<M) {
                                vin[k+UP2] = std::abs(dataMatrixCh1[m+k][n])/NCh1Lin;
                            } else {
                                vin[k+UP2] = 0;
                            }
                        } else {
                            if(n+k>=0 && n+k<Nfull) {
                                rin[k+UP2] = std::abs(dataMatrixCh2[m][n+k])/NCh2Lin;
                            } else {
                                rin[k+UP2] = 0;
                            }
                            if(m+k>=0 && m+k<M) {
                                vin[k+UP2] = std::abs(dataMatrixCh2[m+k][n])/NCh2Lin;
                            } else {
                                vin[k+UP2] = 0;
                            }
                        }
                    }

                    // Range and velocity increased accuracy via FFT interpolation
                    if(upsampleFT(rout, rin, UpsampleFactor)<0) {
                        std::cerr << "Error in upsample of range" << std::endl;
                        return;
                    }
                    std::vector<PRECISION>::iterator rresult = std::max_element(rout.begin(), rout.end());
                    long int rangeMax = std::distance(rout.begin(), rresult);
                    double rDiff = (((double) rangeMax / UpsampleFactor) - UP2)*dr;

                    if(upsampleFT(vout, vin, UpsampleFactor)<0) {
                        std::cerr << "Error in upsample of velocity" << std::endl;
                        return;
                    }
                    std::vector<PRECISION>::iterator vresult = std::max_element(vout.begin(), vout.end());
                    long int velocMax = std::distance(vout.begin(), vresult);
                    double vDiff = -(((double) velocMax / UpsampleFactor) - UP2)*dv;

                    // Storing detection
                    win->Detections->azimuth.push_back(det.azimuth[SNRMaxInd]);
                    win->Detections->range.push_back(det.range[SNRMaxInd]+rDiff);
                    win->Detections->velocity.push_back(det.velocity[SNRMaxInd]+vDiff);
                    win->Detections->snr.push_back(det.snr[SNRMaxInd]);
                    win->Detections->diff_phase.push_back(det.diff_phase[SNRMaxInd]);
                    win->Detections->rindices.push_back(det.rindices[SNRMaxInd]);
                    win->Detections->vindices.push_back(det.vindices[SNRMaxInd]);

                    std::vector<long int> detIndices;
                    // Finding detections "close" to the main peak
                    for(long int i=0; i<det.snr.size(); i++) {
                        // Check if distance to largest detection is smaller than waveform (in range and doppler)
                        if(std::sqrt(std::pow(det.range[i]-det.range[SNRMaxInd],2))<=rWFSize && std::sqrt(std::pow(det.velocity[i]-det.velocity[SNRMaxInd],2))<=vWFSize) {
                            detIndices.push_back(i);
                        }
                    }
                    det.remove_detections(detIndices);
                }

                boost::chrono::system_clock::time_point det_time_stop = boost::chrono::system_clock::now();
                boost::chrono::duration<double> det_time_dur = det_time_stop - det_time_start;

                std::cout << " - processing lasted " << det_time_dur.count() << " seconds" << std::endl;

                // Updating PPI, if activated
                win->emitSignalUpdatePPI(f0);
            }

            { // Tracking
                double phaseDiff = 0;
                boost::chrono::duration<double> t_track = boost::chrono::system_clock::now() - win->start;
                for(long int i=0; i<win->trackers.size(); i++) {
                    // Track motion update
                    win->trackers[i]->motionUpdate(t_track.count(), true);

                    if(win->Detections->range.size()>0){
                        // Get predictions
                        vector<double> cov = win->trackers[i]->covariance();
                        double r,a,v,s;
                        win->trackers[i]->trackPrediction(r,a,v,s);

                        int closestDetection=0;
                        double drMin = absval(win->Detections->range[closestDetection]-r);
                        double dvMin = absval(win->Detections->velocity[closestDetection] - v);
                        double deltar,deltav;
                        for(long int k=1; k<win->Detections->range.size(); k++) {
                            deltar = absval(win->Detections->range[k]-r);
                            deltav = absval(win->Detections->velocity[k] - v);

                            if(deltar<drMin && deltav<dvMin) {
                                closestDetection = k;
                                drMin=deltar;
                                dvMin=deltav;
                            }
                        }

                        phaseDiff = win->Detections->diff_phase[closestDetection];

                        vector<double> trackDistanceInSigmas = win->trackers[i]->distanceInSigmas(win->Detections->range[closestDetection],win->Detections->velocity[closestDetection],
                                                                                                  win->Detections->snr[closestDetection],dr,dv);

                        if(trackDistanceInSigmas(0)<4 && trackDistanceInSigmas(1)<4) {
                            std::cout << "Track #" << i << " is associated with detection " << closestDetection << ", performing informationUpdate" << std::endl;
                            win->trackers[i]->informationUpdate(win->Detections->range[closestDetection], win->Detections->azimuth[closestDetection],
                                                                win->Detections->velocity[closestDetection], win->Detections->snr[closestDetection], dr, dv);
                            win->trackers[i]->trackEstimate(r,a,v,s);
                            std::cout << "Track #" << i << " estimate: " << r << " m, " << a << " rad, " << v << " m/s, " << s << " dB" << std::endl;
                            cov = win->trackers[i]->covariance();
                            std::cout << "Track #" << i << " has cov: " << std::sqrt(cov(0)) << "m, " << std::sqrt(cov(1)) << "m/s, " << std::sqrt(cov(2)) << "rad, " << std::sqrt(cov(3)) << "dB" << std::endl;

                            win->trackers[i]->closest_detection.range.push_back(win->Detections->range[closestDetection]);
                            win->trackers[i]->closest_detection.velocity.push_back(win->Detections->velocity[closestDetection]);
                            win->trackers[i]->closest_detection.azimuth.push_back(win->Detections->azimuth[closestDetection]);
                            win->trackers[i]->closest_detection.snr.push_back(win->Detections->snr[closestDetection]);
                            win->trackers[i]->closest_detection.diff_phase.push_back(win->Detections->diff_phase[closestDetection]);
                        }else {
                            win->trackers[i]->closest_detection.range.push_back(std::nan(""));
                            win->trackers[i]->closest_detection.velocity.push_back(std::nan(""));
                            win->trackers[i]->closest_detection.azimuth.push_back(std::nan(""));
                            win->trackers[i]->closest_detection.snr.push_back(std::nan(""));
                            win->trackers[i]->closest_detection.diff_phase.push_back(std::nan(""));
                        }
                    }
                }
                for(long int i=(int)win->trackers.size()-1; i>=0; i--) {
                    if(win->trackers[i]->updatesWithoutDetections>4) {
                        delete win->trackers[i];
                        win->trackers.erase(win->trackers.begin()+i);
                    }
                }
                if(win->trackers.size()>0) {
                    vector<double> cov = win->trackers[0]->covariance();
                    win->rcovV.push_back(std::sqrt(cov(0)));
                    win->vcovV.push_back(std::sqrt(cov(2)));
                    win->phDiffV.push_back(phaseDiff);
                } else {
                    win->rcovV.push_back(0);
                    win->vcovV.push_back(0);
                    win->phDiffV.push_back(0);
                }
            }


            {// Storing data for display
                std::cout << "Storing data for diplay";
                boost::chrono::system_clock::time_point dopp_start = boost::chrono::system_clock::now();

                double maxVel = win->ui->maxVelocity->text().toDouble();
                double maxRange = win->ui->maxRange->text().toDouble();

                long int minM = std::min((long int) (maxVel/dv), (long int) (M/2));
                maxVel = minM*dv;
                long int maxM = M - minM;

                long int maxN = std::min((size_t) (maxRange/dr), N);

                long int numRows = minM*2-1;

                for(size_t i=0; i<minM; i++){
                    size_t j = i + minM;
                    win->dataMatrixCh1[j] = std::vector<std::complex<PRECISION>>(maxN);
                    win->dataMatrixCh2[j] = std::vector<std::complex<PRECISION>>(maxN);
                    memcpy(&(win->dataMatrixCh1[j][0]), &(dataMatrixCh1[i][0]), sizeof(std::complex<PRECISION>)*maxN);
                    memcpy(&(win->dataMatrixCh2[j][0]), &(dataMatrixCh2[i][0]), sizeof(std::complex<PRECISION>)*maxN);
                }
                for(size_t i=maxM; i<M; i++){
                    size_t j = i - maxM;
                    win->dataMatrixCh1[j] = std::vector<std::complex<PRECISION>>(maxN);
                    win->dataMatrixCh2[j] = std::vector<std::complex<PRECISION>>(maxN);
                    memcpy(&(win->dataMatrixCh1[j][0]), &(dataMatrixCh1[i][0]), sizeof(std::complex<PRECISION>)*maxN);
                    memcpy(&(win->dataMatrixCh2[j][0]), &(dataMatrixCh2[i][0]), sizeof(std::complex<PRECISION>)*maxN);
                }
                win->dimensions[0] = numRows;
                win->dimensions[1] = maxN;

                boost::chrono::system_clock::time_point dopp_stop = boost::chrono::system_clock::now();
                boost::chrono::duration<double> dopp_dur = dopp_stop - dopp_start;

                std::cout << " - processing lasted " << dopp_dur.count() << " seconds" << std::endl;
            }


            {// Call signal to gui to plot
                win->plotFinished=false;
                win->emitSignalcallPlotFunction();
            }

            // Saving data to file
            if(saveData) {
                std::stringstream fname;

                std::time_t sec, nanosec;
                std::tm *tm;
                if(win->simulation) {
                    sec = boost::chrono::system_clock::to_time_t(last);
                    nanosec = boost::chrono::duration_cast<boost::chrono::nanoseconds>(last-boost::chrono::system_clock::from_time_t(sec)).count();
                } else {
                    sec = (std::time_t) rxTimeGPS;
                    nanosec = ((double) (rxTimeGPS-sec))*1e9;
                }

                tm = std::localtime(&sec);

                char strTime[30];
                strftime(strTime, 30, "%Y%b%d-%H%M%S", tm);
                std::cout << strTime << " = " << rxTimeGPS << std::endl;
                fname << win->folder_name << strTime << "." << (int)(nanosec*1e-3)<< ".dat";

                std::cout << "Writing data to: " << fname.str() << std::endl;

                struct fileheader fh;
                fh.sec = sec;
                fh.nanosec = nanosec;
                fh.PRI = PRI;
                fh.f0 = f0;
                fh.fs = rx_rate;
                fh.bandwidth = BW;
                fh.pulse_length = pulse_length;
                fh.M = M;
                fh.N = N;

                std::ofstream filehandle(fname.str(), std::ios::out | std::ios::binary);
                if(filehandle.is_open()){ // Writing data
                    filehandle.write((char*) &fh, sizeof(struct fileheader));
                    if(win->waveform==NULL) {
                        std::cerr << "waveform variable deleted!!" << std::endl;
                    }
                    filehandle.write((char*) &(win->waveform->front()), sizeof(std::complex<PRECISION>)*N);
                    for(long int i=0; i<M; i++) {
                        filehandle.write((char*) &(saveDataCh1[i]->front()), sizeof(std::complex<PRECISION>)*N);
                        delete saveDataCh1[i];
                    }

                    for(long int i=0; i<M; i++) {
                        filehandle.write((char*) &(saveDataCh2[i]->front()), sizeof(std::complex<PRECISION>)*N);
                        delete saveDataCh2[i];
                    }
                    size_t Ndetections = win->Detections->range.size();
                    filehandle.write((char*) &Ndetections, sizeof(size_t));
                    for(long int i=0; i<Ndetections; i++) {
                        filehandle.write((char*) &(win->Detections->range[i]), sizeof(double));
                        filehandle.write((char*) &(win->Detections->azimuth[i]), sizeof(double));
                        filehandle.write((char*) &(win->Detections->velocity[i]), sizeof(double));
                        filehandle.write((char*) &(win->Detections->snr[i]), sizeof(double));
                    }

                    size_t Ntracks = win->trackers.size();
                    filehandle.write((char*) &Ntracks, sizeof(size_t));
                    for(long int i=0; i<Ntracks; i++) {
                        double r,a,v,s;
                        if(win->trackers[i]->updatesWithoutDetections>0)
                            win->trackers[i]->trackEstimate(r,a,v,s);
                        else
                            win->trackers[i]->trackPrediction(r,a,v,s);

                        filehandle.write((char*) &r, sizeof(double));
                        filehandle.write((char*) &a, sizeof(double));
                        filehandle.write((char*) &v, sizeof(double));
                        filehandle.write((char*) &s, sizeof(double));
                        vector<double> cov = win->trackers[i]->covariance();
                        filehandle.write((char*) &(cov[0]), sizeof(double)*4);
                    }
                } else {
                    std::cerr << "Could not write files to folder: " << fname.str() << std::endl;
                }
                filehandle.close();
            }
        }

        // If simulation, sleep to emulate an update rate
        if(win->simulation) {
            now = boost::chrono::system_clock::now();
            sec = now - last;
            double sleep_time = 0.5-sec.count();
            if(sleep_time>0) {
                boost::this_thread::sleep(boost::posix_time::milliseconds((int)(sleep_time*1e3)));
                std::cout << "Simulation: sleeping to emulate longer update rate, " <<
                             sleep_time << " s" << std::endl;
            }
        }

        while((!win->plotFinished) && win->running) {
            boost::this_thread::sleep(boost::posix_time::milliseconds(10));
        }

        delete rx_buff_ch1;
        rx_buff_ch1 = NULL;
        delete rx_buff_ch2;
        rx_buff_ch2 = NULL;
        delete tx_buff;
        tx_buff = NULL;
        delete win->waveform;
        win->waveform=NULL;
        delete win->waveformSpectrum;
        win->waveformSpectrum=NULL;

        CPINum++;
    }
    if(win->waveform!=NULL) {
        delete win->waveform;
        win->waveform=NULL;
    }
    if(win->waveformSpectrum!=NULL) {
        delete win->waveformSpectrum;
        win->waveformSpectrum=NULL;
    }
}

void MainWindow::usrpConnect()
{
    usrp = uhd::usrp::multi_usrp::make(ui->deviceArgs->currentText().toStdString());

    // Specify subdevices on USRP (Important to set up correc transmitters and receivers), is currently set up for UBX-160 on transmit at RF0, and twinRX for 2 receivers at RF1
    std::cout << "Setting subdevs" << std::endl;
    usrp->set_rx_subdev_spec(uhd::usrp::subdev_spec_t("B:0 B:1"));
    usrp->set_tx_subdev_spec(uhd::usrp::subdev_spec_t("A:0"));
    std::cout << "Setting antennas" << std::endl;
    usrp->set_rx_antenna("RX1",0);
    usrp->set_rx_antenna("RX2",1);
    usrp->set_tx_antenna("TX/RX",0);

    // Set clock source
    usrp->set_clock_source("gpsdo");
    usrp->set_time_source("gpsdo");

    // making sure that ref is locked
    std::cout << "Locking to ref" << std::endl;
    bool locked=false;
    for(int i=0; i<1; i++) {
        if(usrp->get_mboard_sensor("ref_locked").to_bool()) {
            locked=true;
            break;
        }
        else
            boost::this_thread::sleep(boost::posix_time::seconds(1));
    }
    std::cout << usrp->get_mboard_sensor("ref_locked").to_pp_string() << std::endl;
    if(!locked) {
        throw std::runtime_error("Ref not locked");
    }

    // making sure that gps is locked
    std::cout << "Locking to GPS" << std::endl;
    locked=false;
    for(int i=0; i<1; i++) {
        if(usrp->get_mboard_sensor("gps_locked").to_bool()) {
            locked=true;
            break;
        }
        else
            boost::this_thread::sleep(boost::posix_time::seconds(1));
    }
    if(!locked) {
        //throw std::runtime_error("GPS not locked");
        std::cout << "Warning: GPS not locked" << std::endl;
    }
    std::cout << usrp->get_mboard_sensor("gps_locked").to_pp_string() << std::endl;

    // Set GPS device time
    uhd::time_spec_t gps_time = uhd::time_spec_t(time_t(usrp->get_mboard_sensor("gps_time", 0).to_int()));
    usrp->set_time_next_pps(gps_time+1.0, 0);

    // Wait for it to apply
    boost::this_thread::sleep(boost::posix_time::seconds(2));

    // Printf mboard information
    std::cout << usrp->get_pp_string() << std::endl;

    // Setting channel 1 and 2 as companion (twinRX specific)
    usrp->set_rx_lo_source("internal", uhd::usrp::multi_usrp::ALL_LOS, 0);
    usrp->set_rx_lo_source("companion", uhd::usrp::multi_usrp::ALL_LOS, 1);

    // Checking lo config on twinRX
    for(int k=0; k<2; k++) {
        std::cout << "Current LO source Channel " << k << ": " << usrp->get_rx_lo_source(uhd::usrp::multi_usrp::ALL_LOS, k) << std::endl;
    }

    // Sleep 1 second to allow all settings to be done
    boost::this_thread::sleep(boost::posix_time::seconds(1));
}

void MainWindow::usrpUpdateTextFields() {
    double rx_gain, tx_gain, rx_freq, rx_rate, tx_rate;

    rx_gain = usrp->get_rx_gain(0);
    tx_gain = usrp->get_tx_gain(0);
    rx_rate = usrp->get_rx_rate(0);
    tx_rate = usrp->get_tx_rate(0);
    rx_freq = usrp->get_rx_freq(0);

    ui->GainRX->setText(QString::number(rx_gain));
    ui->GainTX->setText(QString::number(tx_gain));
    ui->rxRate->setText(QString::number(rx_rate));
    ui->txRate->setText(QString::number(tx_rate));
    ui->carrierFrequency->setText(QString::number(rx_freq));
}

void MainWindow::usrpSetParametersFromFields() {
    double rx_gain, tx_gain, rx_freq, rx_rate, tx_rate;

    rx_gain = ui->GainRX->text().toDouble();
    tx_gain = ui->GainTX->text().toDouble();
    rx_freq = ui->carrierFrequency->text().toDouble();
    rx_rate = ui->rxRate->text().toDouble();
    tx_rate = ui->txRate->text().toDouble();

    // Setting rate
    usrp->set_rx_rate(rx_rate, 0);
    usrp->set_rx_rate(rx_rate, 1);
    usrp->set_tx_rate(tx_rate, 0);

    // Setting gain
    usrp->set_rx_gain(rx_gain, 0);
    usrp->set_rx_gain(rx_gain, 1);
    usrp->set_tx_gain(tx_gain, 0);

    // SEtting freq
    usrp->set_rx_freq(uhd::tune_request_t(rx_freq), 0);
    usrp->set_rx_freq(uhd::tune_request_t(rx_freq), 1);
    usrp->set_tx_freq(uhd::tune_request_t(rx_freq, 80e6), 0);
    { // Timed set freq
        uhd::time_spec_t cmd_time = usrp->get_time_now() + uhd::time_spec_t(0.1);
        usrp->set_command_time(cmd_time);
        usrp->set_rx_freq(uhd::tune_request_t(rx_freq), 0);
        usrp->set_rx_freq(uhd::tune_request_t(rx_freq), 1);
        usrp->set_tx_freq(uhd::tune_request_t(rx_freq, 80e6), 0);
        usrp->clear_command_time();
    }
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    running=false;
    ppi=false;
    dataMatrixCh1 = std::vector<std::vector<std::complex<PRECISION>>>(MAX_NUM_PULSES);
    dataMatrixCh2 = std::vector<std::vector<std::complex<PRECISION>>>(MAX_NUM_PULSES);
    dimensions = std::vector<size_t>(2,0);

    Detections = new Detection;

    ppiWindow = new ppiDialog(this);
    cameraWindow = new cameraDialog(this);

    connect(ui->qwtPlot, SIGNAL(selectedPoint(QPointF)),
            this, SLOT(qwtPlotPointSelected(QPointF)));

    connect(this, SIGNAL(callPlotFunction()), this, SLOT(plotFunction()));
    connect(this, SIGNAL(callUpdatePPI(Detection*, double)), ppiWindow, SLOT(updatePPI(Detection*, double)));

    this->statusBar()->show();

    start = boost::chrono::system_clock::now();

    w_curve = new QwtPlotCurve();
    w_curve->setStyle(QwtPlotCurve::Lines);
    w_curve->attach(ui->qwtWaveform);

    ws_curve = new QwtPlotCurve();
    ws_curve->setStyle(QwtPlotCurve::Lines);
    ws_curve->attach(ui->qwtWaveformSpectrum);

    prf_curve = new QwtPlotCurve();
    prf_curve->setStyle(QwtPlotCurve::Lines);
    prf_curve->attach(ui->qwtPRF);

    bw_curve = new QwtPlotCurve();
    bw_curve->setStyle(QwtPlotCurve::Lines);
    bw_curve->attach(ui->qwtBandwidth);

    pl_curve = new QwtPlotCurve();
    pl_curve->setStyle(QwtPlotCurve::Lines);
    pl_curve->attach(ui->qwtPulselength);

    numpulses_curve = new QwtPlotCurve();
    numpulses_curve->setStyle(QwtPlotCurve::Lines);
    numpulses_curve->attach(ui->qwtNumPulses);

    plotFinished=true;

    ui->txtFileName->setText(QString(getenv("HOME")) + QString("/DATA/Waveforms/"));

    waveform = NULL;
    waveformSpectrum = NULL;
}

MainWindow::~MainWindow()
{
    delete ppiWindow;
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    this->statusBar()->showMessage("Connecting to USRP...");
    if(ui->deviceArgs->currentText() == QString("Simulation")) {
        simulation = true;
        targets.push_back(new SimulatedTarget(300, 600, -20, 0, 0.1));
    } else {
        simulation = false;

        usrpConnect();
        usrpSetParametersFromFields();
        usrpUpdateTextFields();
    }
    if(usrp==NULL && !simulation) {
        this->statusBar()->showMessage("Could not connect");
    } else {
        this->statusBar()->showMessage("Connected!!");
    }
}

void MainWindow::on_pushButton_3_clicked()
{
    if(!simulation) {
        usrpSetParametersFromFields();
        usrpUpdateTextFields();
    }
}

void MainWindow::on_pushButton_2_clicked()
{
    running = !running;

    ui->pushButton_2->setChecked(running);

    if(running) {
        firstPlot=true;
        for(int i=0; i<trackers.size(); i++)
            delete trackers[i];
        trackers.clear();
        pulse_doppler_thread.create_thread(boost::bind(&pulse_doppler_worker, this));
    } else {
        pulse_doppler_thread.join_all();
    }
}

void MainWindow::plotFunction()
{
    if(running) {
        double rx_rate = ui->rxRate->text().toDouble();

        // Waveform plot
        std::vector<double> t(waveformDim), wf(waveformDim);
        for(long int k=0; k<waveformDim; k++) {
            wf[k] = std::real((*waveform)[k]);
            t[k] = k/rx_rate;
        }
        w_curve->setSamples(t.data(), wf.data(), (int)waveformDim);
        ui->qwtWaveform->replot();

        // Waveform spectrum plot
        std::vector<double> f(waveformSpectrumDim), wfs(waveformSpectrumDim);
        for(long int k=0; k<waveformSpectrumDim; k++) {
            if(k<std::floor(waveformSpectrumDim/2)) {
                wfs[k] = 20*std::log10(std::abs((*waveformSpectrum)[k+std::floor(waveformSpectrumDim/2)]));
            } else {
                wfs[k] = 20*std::log10(std::abs((*waveformSpectrum)[k-std::floor(waveformSpectrumDim/2)]));
            }
            f[k] = 1e-6*(k-std::floor(waveformSpectrumDim/2))*rx_rate/waveformSpectrumDim;
        }
        ws_curve->setSamples(f.data(), wfs.data(), (int)waveformSpectrumDim);
        ui->qwtWaveformSpectrum->replot();

        // PRF plot
        prf_curve->setSamples(tV.data(), prfV.data(), (int)tV.size());
        ui->qwtPRF->replot();

        // Bandwidth plot
        bw_curve->setSamples(tV.data(), bwV.data(), (int)tV.size());
        ui->qwtBandwidth->replot();

        // Pulse length plot
        pl_curve->setSamples(tV.data(), plV.data(), (int)tV.size());
        ui->qwtPulselength->replot();

        // Pulse length plot
        numpulses_curve->setSamples(tV.data(), npV.data(), (int)tV.size());
        ui->qwtNumPulses->replot();

        // Detection plots
        ui->qwtPlot->d_curve->setSamples(Detections->range.data(), Detections->velocity.data(), (int)Detections->range.size());

        // Track plots
        std::vector<double> R,V;
        double r,a,v,s;
        for(long int i=0; i<trackers.size(); i++) {
            trackers[i]->trackEstimate(r,a,v,s);
            R.push_back(r);
            V.push_back(v);
        }
        ui->qwtPlot->t_curve->setSamples(R.data(), V.data(), (int) R.size());

        // Range-Doppler plot
        double maxVel = ui->maxVelocity->text().toDouble();
        double maxRange = ui->maxRange->text().toDouble();

        QVector<double> *chRDMap = new QVector<double>(dimensions[0]*dimensions[1]);

        if( ui->channelPlot->currentIndex() == PLOT_CHANNEL_1 ) {
            for(long int i=0; i<dimensions[0]; i++) {
                for(long int k=0; k<dimensions[1]; k++) {
                    // Copying data to vector backwards to get v=-fd/lambda relationship
                    (*chRDMap)[i*dimensions[1]+k] = 20*std::log10(std::abs(dataMatrixCh1[dimensions[0]-i][k])) - NoiseFloorCh1;
                }
            }
        } else {
            for(long int i=0; i<dimensions[0]; i++) {
                for(long int k=0; k<dimensions[1]; k++) {
                    (*chRDMap)[i*dimensions[1]+k] = 20*std::log10(std::abs(dataMatrixCh2[dimensions[0]-i][k])) - NoiseFloorCh1;
                }
            }
        }
        QwtMatrixRasterData *ch1Raster = new QwtMatrixRasterData();

        ch1Raster->setValueMatrix((*chRDMap), dimensions[1]);

        ch1Raster->setInterval( Qt::YAxis, QwtInterval( -maxVel, maxVel, QwtInterval::ExcludeMaximum ) );
        ch1Raster->setInterval( Qt::XAxis, QwtInterval( 0, maxRange, QwtInterval::ExcludeMaximum ) );
        ch1Raster->setInterval( Qt::ZAxis, QwtInterval( ui->signalMin->text().toDouble(), ui->signalMax->text().toDouble() ) );

        ui->qwtPlot->d_spectrogram->setData(ch1Raster);

        QwtInterval zInterval = ui->qwtPlot->d_spectrogram->data()->interval( Qt::ZAxis );
        // A color bar on the right axis
        QwtScaleWidget *rightAxis = ui->qwtPlot->axisWidget( QwtPlot::yRight );
        rightAxis->setColorBarEnabled( true );
        rightAxis->setColorBarWidth( 40 );
        rightAxis->setColorMap( zInterval, new ColorMap() );

        ui->qwtPlot->setAxisScale( QwtPlot::yRight, zInterval.minValue(), zInterval.maxValue() );
        ui->qwtPlot->enableAxis( QwtPlot::yRight );

        ui->qwtPlot->plotLayout()->setAlignCanvasToScales( true );

        ui->qwtPlot->setAxisScale( QwtPlot::yLeft, -maxVel, maxVel);
        ui->qwtPlot->setAxisMaxMinor( QwtPlot::yLeft, 0 );
        ui->qwtPlot->setAxisScale( QwtPlot::xBottom, 0, maxRange);
        ui->qwtPlot->setAxisMaxMinor( QwtPlot::xBottom, 0 );

        // Replotting
        ui->qwtPlot->replot();

        delete chRDMap;
    }
    plotFinished = true;
}

void MainWindow::qwtPlotPointSelected(QPointF point) {
    if(Detections->range.size()<1) {
        qDebug() << "No detections, cannot initiate track";
        return;
    }

    double range = point.x();
    double velocity = point.y();

    long int closestIndex=0;
    for(long int i=1; i<Detections->range.size(); i++) {
        if(std::abs(Detections->range[i]-range) < std::abs(Detections->range[closestIndex]-range) &&
                std::abs(Detections->velocity[i]-velocity) < std::abs(Detections->velocity[closestIndex]-velocity)) {
            closestIndex=i;
        }
    }

    double drMin = 500;
    double dvMin = 10;
    if((std::abs(Detections->range[closestIndex]-range)>drMin || std::abs(Detections->velocity[closestIndex]-velocity)>dvMin)) {
        qDebug() << "Centerclick too far from detection, disregarding...";
        return;
    }

    qDebug() << "Starting track at: " << Detections->range[closestIndex] << " m, " <<
                Detections->velocity[closestIndex] << " m/s, " << Detections->snr[closestIndex] << " dB";


    boost::chrono::system_clock::time_point now = boost::chrono::system_clock::now();
    boost::chrono::duration<double> sec = now - start;

    double P_0 = 0.5;
    double P_max = 0.03;
    double A_max = 3*9.81;
    double alpha = (double) 1/30;

    Tracker *track = new Tracker(sec.count(), Detections->range[closestIndex], Detections->velocity[closestIndex], Detections->azimuth[closestIndex], Detections->snr[closestIndex],
                                       alpha, P_0, P_max, A_max);
    trackers.push_back(track);
}

void MainWindow::on_saveDataToFile_stateChanged(int arg1)
{
    if(arg1>0) {
        time_t tnow = boost::chrono::system_clock::to_time_t( boost::chrono::system_clock::now() );
        struct tm * now = localtime( & tnow );

        std::ostringstream strs;

        char strTime[30];
        strftime(strTime, 30, "%Y%b%d-%H%M%S", now);
        strs << getenv("HOME") << "/DATA/" << strTime << "/";

        folder_name = strs.str();
        boost::filesystem::path dir(folder_name.c_str());
        boost::filesystem::create_directory(dir);

        if(camera) {
            const std::string videoName = folder_name + std::string(strTime) + "_video.avi";
            cameraWindow->start_video_writer(videoName);
        }
    } else {
        if(camera) {
            cameraWindow->stop_video_writer();
        }
    }
}

void MainWindow::on_actionSync_to_GPS_triggered()
{
    // Check if GPS is locked
    std::cout << "Locking to GPS" << std::endl;
    bool locked=false;
    for(int i=0; i<1; i++) {
        if(usrp->get_mboard_sensor("gps_locked").to_bool()) {
            locked=true;
            break;
        }
        else
            boost::this_thread::sleep(boost::posix_time::seconds(1));
    }
    if(!locked) {
        std::cout << "Warning: GPS not locked" << std::endl;
    }
    std::cout << usrp->get_mboard_sensor("gps_locked").to_pp_string() << std::endl;

    // Set GPS device time
    uhd::time_spec_t gps_time = uhd::time_spec_t(time_t(usrp->get_mboard_sensor("gps_time", 0).to_int()));
    usrp->set_time_next_pps(gps_time+1.0, 0);
}

void MainWindow::on_pushButton_4_clicked()
{
    ppi = !ppi;

    ui->pushButton_4->setChecked(ppi);

    if(ppi) {
        ppiWindow->show();
    } else {
        ppiWindow->hide();
    }
}

void MainWindow::on_numPulses_editingFinished()
{
    long int M = ui->numPulses->text().toLong();
    long int Mmod2 = std::pow(2,std::ceil(std::log2(M)));
    if(M!=Mmod2) {
        ui->numPulses->setText(QString::number(Mmod2));
        std::cout << "NumPulses set is not 2 multiple" << std::endl;
    }
}

void MainWindow::on_checkBox_toggled(bool checked)
{
    camera = checked;
    cameraWindow->m_enable_camera = camera;

    if(checked) {
        // Show camera dialog
        cameraWindow->show();
//        cameraWindow->start_timer();
    } else {
        // Close camera dialog
        cameraWindow->hide();
        cameraWindow->stop_timer();
    }
}

void MainWindow::on_numPulses2_editingFinished()
{
    //ui->numPulses->setText(QString::number(std::pow(2, ui->numPulses2->value())));
}

void MainWindow::on_numPulses2_valueChanged(int arg1)
{
    ui->numPulses->setText(QString::number(std::pow(2, arg1)));
}

void MainWindow::on_actionConfig_triggered(bool checked)
{
    checked ? ui->frConfig->show() : ui->frConfig->hide();
}

