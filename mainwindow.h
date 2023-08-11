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
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTextStream>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <random>
#include <ctime>

#include <uhd/types/tune_request.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>

#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_zoomer.h>
#include <qwt/qwt_symbol.h>
#include <qwt/qwt_legend.h>
#include <qwt/qwt_matrix_raster_data.h>
#include <qwt/qwt_scale_widget.h>

#include <boost/assign.hpp>
#include <boost/thread/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/chrono/chrono.hpp>
#include <boost/filesystem.hpp>
#include <boost/circular_buffer.hpp>

#include "cameradialog.h"
#include "plot.h"
#include "simulatedtarget.h"
#include "tracker.h"
#include "processing.h"
#include "ppidialog.h"
#include "detection.h"
#include "waveform.h"

#define SAVE_TYPE_RANGE_DOPPLER 0
#define SAVE_TYPE_PULSE_COMPRESSION 1
#define SAVE_TYPE_RAW 2

#define PLOT_CHANNEL_1 0
#define PLOT_CHANNEL_2 1

#define MAX_NUM_PULSES 100000
#define MAX_NUM_SAMPLES 500000000

#define MIN_PRI 1e-6
#define MAX_PRI 1e-2
#define MIN_CPI 1e-6
#define MAX_CPI 4

struct fileheader {
    size_t sec;
    size_t nanosec;
    double PRI;
    double pulse_length;
    double f0;
    double fs;
    double bandwidth;
    size_t M;
    size_t N;
};


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    bool running, ppi, plotFinished, camera;
    uhd::usrp::multi_usrp::sptr usrp;
    Ui::MainWindow *ui;
    std::vector<std::vector<std::complex<PRECISION>>> dataMatrixCh1;
    std::vector<std::vector<std::complex<PRECISION>>> dataMatrixCh2;
    double NoiseFloorCh1;
    double NoiseFloorCh2;
    std::vector<std::complex<PRECISION>> *waveformSpectrum, *waveform;
    size_t waveformSpectrumDim, waveformDim;
    std::vector<size_t> dimensions;
    Detection *Detections;
    bool simulation;
    std::vector<SimulatedTarget*> targets;
    std::vector<Tracker*> trackers;
    boost::chrono::system_clock::time_point start;
    std::string folder_name;
    std::vector<double> prfV, bwV, plV, npV, rcovV, vcovV, tV, phDiffV;
    QwtPlotCurve *ws_curve, *w_curve, *prf_curve, *bw_curve, *pl_curve, *numpulses_curve, *rcov_curve, *vcov_curve, *phd_curve;

    ppiDialog *ppiWindow;
    cameraDialog * cameraWindow;

    void usrpSetParametersFromFields();

    void emitSignalcallPlotFunction() {
        emit callPlotFunction();
    }

    void emitSignalUpdatePPI(double f0) {
        if(ppi) {
            emit callUpdatePPI(Detections, f0);
        }
    }


signals:
    void callPlotFunction();
    void callUpdatePPI(Detection *, double);

private slots:
    void on_pushButton_clicked();

    void on_pushButton_3_clicked();

    void on_pushButton_2_clicked();

    void qwtPlotPointSelected(QPointF point);
    void plotFunction();

    void on_saveDataToFile_stateChanged(int arg1);

    void on_actionSync_to_GPS_triggered();

    void on_pushButton_4_clicked();

    void on_numPulses_editingFinished();

    void on_checkBox_toggled(bool checked);

    void on_numPulses2_editingFinished();

    void on_numPulses2_valueChanged(int arg1);
    void on_actionConfig_triggered(bool checked);

private:
    boost::thread_group pulse_doppler_thread,camera_thread;
    bool firstPlot, m_enable_camera;

    void usrpConnect();
    void usrpUpdateTextFields();

};

#endif // MAINWINDOW_H
