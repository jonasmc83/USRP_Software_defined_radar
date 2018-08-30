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
#ifndef PPIDIALOG_H
#define PPIDIALOG_H

#include <QDialog>
#include "detection.h"
#include <cmath>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>

namespace Ui {
class ppiDialog;
}

class ppiDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ppiDialog(QWidget *parent = 0);
    ~ppiDialog();
    double phaseToAzimuth(double phase, double f0);

private:
    Ui::ppiDialog *ui;

    QwtPlotCurve *dets, *hist;
    int historyLen;
    std::vector<std::vector<double>*> xHistory;
    std::vector<std::vector<double>*> yHistory;

private slots:
    void updatePPI(Detection *, double);
    void on_northRange_editingFinished();
    void on_eastRange_editingFinished();
    void on_historyLength_valueChanged(int arg1);
};

#endif // PPIDIALOG_H
