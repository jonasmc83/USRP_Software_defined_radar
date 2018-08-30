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
#ifndef PLOT_H
#define PLOT_H

//#define QWT_DLL

#include <qwt/qwt_plot.h>
#include <qwt/qwt_plot_spectrogram.h>
#include <qwt/qwt_color_map.h>
#include <qwt/qwt_plot_spectrogram.h>
#include <qwt/qwt_plot_layout.h>
#include <qwt/qwt_matrix_raster_data.h>
#include <qwt/qwt_scale_widget.h>
#include <qwt/qwt_plot_magnifier.h>
#include <qwt/qwt_plot_panner.h>
#include <qwt/qwt_plot_renderer.h>
#include <qwt/qwt_plot_grid.h>
#include <qwt/qwt_plot_canvas.h>
#include <qwt/qwt_plot_picker.h>
#include <qwt/qwt_picker_machine.h>
#include <qwt/qwt_plot_curve.h>
#include <qwt/qwt_symbol.h>

class ColorMap: public QwtLinearColorMap
{
public:
    ColorMap():
        QwtLinearColorMap( Qt::darkBlue, Qt::darkRed )
    {
        addColorStop( 0.2, Qt::blue );
        addColorStop( 0.4, Qt::cyan );
        addColorStop( 0.6, Qt::yellow );
        addColorStop( 0.8, Qt::red );
    }
};

class Plot: public QwtPlot
{
    Q_OBJECT

public:
    Plot( QWidget * = NULL );
    QwtPlotSpectrogram *d_spectrogram;
    QwtPlotCurve *d_curve;
    QwtPlotCurve *t_curve;

public Q_SLOTS:
    void exportPlot();
    void setResampleMode( int );
    void pointSelected(QPointF);

signals:
    void selectedPoint( QPointF );
};

#endif // PLOT_H
