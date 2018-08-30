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
#ifndef CAMERADIALOG_H
#define CAMERADIALOG_H

#include <QDialog>
#include <QTimer>
#include <QPainter>
#include <QThread>

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <iostream>

namespace Ui {
class cameraDialog;
}

class cameraDialog : public QDialog
{
    Q_OBJECT

public:
    explicit cameraDialog(QWidget *parent = 0);
    ~cameraDialog();

    bool m_enable_camera;

    void start_timer();
    void stop_timer();
    void start_video_writer(const std::string fname);
    void stop_video_writer();

private:
    Ui::cameraDialog *ui;

    QImage m_i;

    cv::VideoCapture capture;
    cv::VideoWriter writer;
    cv::Mat frame;

    bool writing;
    int fps;

    QTimer* m_timer;
    QThread *m_thread;

    void paintEvent(QPaintEvent* e);
    QImage Mat2QImage(cv::Mat &src);

public slots:
    void queryFrame();
private slots:
    void on_fpsEdit_editingFinished();
    void on_pushButton_clicked();
};

#endif // CAMERADIALOG_H
