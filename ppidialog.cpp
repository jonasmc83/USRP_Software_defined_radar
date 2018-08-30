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
#include "ppidialog.h"
#include "ui_ppidialog.h"

ppiDialog::ppiDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ppiDialog)
{
    ui->setupUi(this);

    ui->qwtPPI->setAxisScale(QwtPlot::xBottom, -(ui->eastRange->text().toDouble()), ui->eastRange->text().toDouble());
    ui->qwtPPI->setAxisScale(QwtPlot::yLeft, -(ui->northRange->text().toDouble()), ui->northRange->text().toDouble());

    dets = new QwtPlotCurve();
    dets->setStyle(QwtPlotCurve::Dots);
    dets->setSymbol(new QwtSymbol(QwtSymbol::Diamond, QColor(Qt::black), QColor(Qt::red), QSize(7,7)));
    dets->attach(ui->qwtPPI);

    hist = new QwtPlotCurve();
    hist->setStyle(QwtPlotCurve::Dots);
    hist->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, QColor(Qt::black), QColor(Qt::blue), QSize(5,5)));
    hist->attach(ui->qwtPPI);

    xHistory = std::vector<std::vector<double>*>(0);
    yHistory = std::vector<std::vector<double>*>(0);

    historyLen = ui->historyLength->value();
}

ppiDialog::~ppiDialog()
{
    for(int k = 0; k<xHistory.size(); k++) {
        delete xHistory[k];
        delete yHistory[k];
    }
    delete ui;
}

void ppiDialog::updatePPI(Detection *det, double f0) {
    long int M = det->range.size();
    std::vector<double> histX(0), histY(0);
    std::vector<double> *x = new std::vector<double>(M);
    std::vector<double> *y = new std::vector<double>(M);

    for(long int k=0; k<M; k++) {
        double azRad = (M_PI/180)*phaseToAzimuth(det->diff_phase[k], f0);
        (*x)[k] = det->range[k] * std::sin(azRad);
        (*y)[k] = det->range[k] * std::cos(azRad);
    }
    dets->setSamples(x->data(), y->data(), M);

    for(int k=0; k<xHistory.size(); k++) {
        histX.insert(histX.end(), (xHistory[k])->begin(), (xHistory[k])->end());
        histY.insert(histY.end(), (yHistory[k])->begin(), (yHistory[k])->end());
    }
    hist->setSamples(histX.data(), histY.data(), histX.size());

    xHistory.insert(xHistory.begin(), x);
    yHistory.insert(yHistory.begin(), y);

    while(xHistory.size()>historyLen) {
        delete xHistory[xHistory.size()-1];
        delete yHistory[yHistory.size()-1];
        xHistory.pop_back();
        yHistory.pop_back();
    }

    ui->qwtPPI->replot();
}

double ppiDialog::phaseToAzimuth(double phase, double f0) {
    double corr_phase = phase - (M_PI/180)*(ui->receiverOffset->text().toDouble());
    double lambda = 3e8 / f0;
    return (180/M_PI) * std::asin( (lambda * corr_phase) / (2 * M_PI * ui->elementSpacing->text().toDouble()) ) + ui->boresightBearing->text().toDouble();
}

void ppiDialog::on_northRange_editingFinished()
{
    ui->qwtPPI->setAxisScale(QwtPlot::yLeft, -(ui->northRange->text().toDouble()), ui->northRange->text().toDouble());
    ui->qwtPPI->replot();
}

void ppiDialog::on_eastRange_editingFinished()
{
    ui->qwtPPI->setAxisScale(QwtPlot::xBottom, -(ui->eastRange->text().toDouble()), ui->eastRange->text().toDouble());
    ui->qwtPPI->replot();
}


void ppiDialog::on_historyLength_valueChanged(int arg1)
{
    historyLen = arg1;
}
