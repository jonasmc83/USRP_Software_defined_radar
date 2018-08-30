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
#include "plot.h"

Plot::Plot( QWidget *parent ):
    QwtPlot( parent )
{
    QwtPlotCanvas *canvas = new QwtPlotCanvas();
    canvas->setBorderRadius( 10 );
    setCanvas( canvas );

    d_spectrogram = new QwtPlotSpectrogram();
    d_spectrogram->setRenderThreadCount( 0 ); // use system specific thread count
    d_spectrogram->setColorMap( new ColorMap() );
    d_spectrogram->attach( this );

    d_curve = new QwtPlotCurve();
    d_curve->setStyle(QwtPlotCurve::NoCurve);
    d_curve->setSymbol(new QwtSymbol(QwtSymbol::Diamond, Qt::white, QPen(Qt::black, 1), QSize(10,10)));
    d_curve->attach(this);

    t_curve = new QwtPlotCurve();
    t_curve->setStyle(QwtPlotCurve::NoCurve);
    t_curve->setSymbol(new QwtSymbol(QwtSymbol::RTriangle, Qt::green, QPen(Qt::black, 1), QSize(7,7)));
    t_curve->attach(this);

    QwtPlotMagnifier *magnifier = new QwtPlotMagnifier( canvas );
    magnifier->setAxisEnabled( QwtPlot::yRight, false );

    QwtPlotPanner *panner = new QwtPlotPanner( canvas );
    panner->setAxisEnabled( QwtPlot::yRight, false );

    QwtPlotPicker* clickPicker = new QwtPlotPicker(this->canvas());
    clickPicker->setStateMachine(new QwtPickerClickPointMachine);
    clickPicker->setMousePattern(QwtPicker::MouseSelect1, Qt::MiddleButton);
    connect(clickPicker, SIGNAL(appended(QPointF)),
            this, SLOT(pointSelected(QPointF)));
}

void Plot::exportPlot()
{
    QwtPlotRenderer renderer;
    renderer.exportTo( this, "rasterview.pdf" );
}

void Plot::setResampleMode( int mode )
{
    QwtMatrixRasterData *data = static_cast<QwtMatrixRasterData *>( d_spectrogram->data() );
    data->setResampleMode(
        static_cast<QwtMatrixRasterData::ResampleMode>( mode ) );

    replot();
}

void Plot::pointSelected(QPointF point) {
    emit selectedPoint(point);
}
