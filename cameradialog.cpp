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
#include "cameradialog.h"
#include "ui_cameradialog.h"

cameraDialog::cameraDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::cameraDialog)
{
    ui->setupUi(this);

    m_thread = new QThread;
    this->moveToThread(m_thread);

    m_timer = new QTimer(this);
    connect(m_timer,SIGNAL(timeout()),this,SLOT(queryFrame()));

    writing = false;
}

cameraDialog::~cameraDialog()
{
    delete m_thread;
    delete m_timer;
    delete ui;
}

void cameraDialog::queryFrame() {
    if (!capture.read(frame)) {
        std::cout << "Could not read frame" << std::endl;
        return;
    }
    if(writing) {
        if(writer.isOpened()) {
            writer.write(frame);
        } else {
            std::cout << "Writer is not open" << std::endl;
        }
    }

    this->update();
}


void cameraDialog::paintEvent(QPaintEvent* e) {
    QPainter painter(this);

    if(frame.empty())
        return;

    m_i = Mat2QImage(frame);

    painter.drawImage(QPoint(ui->frame->x(),ui->frame->y()),m_i);
}

QImage cameraDialog::Mat2QImage(cv::Mat &src)
{
    cv::Mat temp(src.cols,src.rows,src.type());
    cvtColor(src, temp, cv::COLOR_BGR2RGB);
    QImage dest= QImage((uchar*) temp.data, temp.cols, temp.rows, temp.step, QImage::Format_RGB888);

    return dest;
}

void cameraDialog::start_timer() {
    capture = cv::VideoCapture(ui->lineURL->text().toStdString());

    if (!capture.isOpened()) {
        std::cout << "Could not open url" << std::endl;
        return;
    }

    if (!capture.read(frame)) {
        std::cout << "Could not read frame" << std::endl;
        return;
    }

    fps = ui->fpsEdit->text().toInt();

    //set timer for (1/fps)*1000ms intervals
    m_timer->start(1000/((double)fps));
}

void cameraDialog::stop_timer() {
    capture.release();

    m_timer->stop();
}

void cameraDialog::start_video_writer(const std::string fname) {
    if(!capture.isOpened()) {
        std::cout << "Could not write video since no camera is opened" << std::endl;
        return;
    }

    int codec = cv::VideoWriter::fourcc('M', 'J', 'P', 'G');
    cv::Size sz = cv::Size(1920,1080);

    std::cout << "Opening video file for writing: " << fname << std::endl;
    writer = cv::VideoWriter(fname, codec, fps,  sz);

    if(!writer.isOpened()) {
        std::cout << "Could not open video file for writing" << std::endl;
    }

    writing = true;
}

void cameraDialog::stop_video_writer() {
    if(writer.isOpened()) {
        writer.release();
        writing = false;
    }
}


void cameraDialog::on_fpsEdit_editingFinished()
{
    fps = ui->fpsEdit->text().toInt();
    m_timer->setInterval(1000/((double)fps));
}

void cameraDialog::on_pushButton_clicked()
{
    stop_timer();
    start_timer();
}
