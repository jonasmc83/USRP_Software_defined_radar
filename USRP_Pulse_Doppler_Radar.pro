#    USRP_Software_defined_radar is a software for real time sampling, processing, display and storing
#    Copyright (C) 2018  Jonas Myhre Christiansen <jonas-myhre.christiansen@ffi.no>
#	
#    This file is part of USRP_Software_defined_radar.
#
#    USRP_Software_defined_radar is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    USRP_Software_defined_radar is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with USRP_Software_defined_radar.  If not, see <https://www.gnu.org/licenses/>.

#-------------------------------------------------
#
# Project created by QtCreator 2017-01-11T11:14:29
#
#-------------------------------------------------

QT += core gui
CONFIG += c++11

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = USRP_Pulse_Doppler_Radar
TEMPLATE = app

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    plot.cpp \
    processing.cpp \
    simulatedtarget.cpp \
    tracker.cpp \
    ppidialog.cpp \
    detection.cpp \
    cameradialog.cpp \
    waveform.cpp

HEADERS += \
    mainwindow.h \
    plot.h \
    processing.h \
    simulatedtarget.h \
    tracker.h \
    ppidialog.h \
    detection.h \
    processing_gpu.h \
    processing_gpu_double.h \
    precision.h \
    cameradialog.h \
    waveform.h

FORMS += \
    mainwindow.ui \
    ppidialog.ui \
    cameradialog.ui

#========== Libraries ==========#
# QWT
CONFIG += qwt
DEFINES += QWT_DLL
LIBS += -lqwt-qt5
INCLUDEPATH += /usr/include/qwt

# UHD
LIBS += -luhd -lboost_filesystem -lboost_system -lboost_thread -lboost_chrono

# OpenCV
INCLUDEPATH += /usr/include/opencv4
OPENCV_DIR = /usr/lib/x86_64-linux-gnu
LIBS += -L$$OPENCV_DIR -lopencv_core -lopencv_highgui -lopencv_videoio -lopencv_imgproc

# OpenMP
#LIBS += -lgomp
QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

# KissFFT
INCLUDEPATH += /usr/include/kissfft

#========== CUDA ==========#
#CUDA_DIR = /usr/local/cuda-8.0  # Path to cuda toolkit install
!isEmpty(CUDA_DIR): {
    INCLUDEPATH += $$CUDA_DIR/include
    QMAKE_LIBDIR += $$CUDA_DIR/lib64
    NVCC_PATH = $$CUDA_DIR/bin/nvcc
} else { # System path
    NVCC_PATH = nvcc
}

# Build options
CUDA_OBJECTS_DIR = $$OUT_PWD/cuda
CUDA_ARCH = sm_52 # Type of CUDA architecture
NVCC_OPTIONS = --use_fast_math

# CUDA libraries
LIBS += -lcuda -lcudart -lcufft

# The following makes sure all path names (which often include spaces) are put between quotation marks
CUDA_INC = $$join(INCLUDEPATH,'" -I"','-I"','"')

# Adding CUDA sources
CUDA_SOURCES += \
    processing_gpu.cu \
    processing_gpu_double.cu

OTHER_FILES += $$CUDA_SOURCES

# Configuration of the Cuda compiler
cuda.input = CUDA_SOURCES
cuda.output = $$CUDA_OBJECTS_DIR/${QMAKE_FILE_BASE}_cuda.o
cuda.commands = $$NVCC_PATH $$NVCC_OPTIONS $$CUDA_INC $$LIBS -m64 -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
CONFIG(debug, debug|release): cuda.commands += -D_DEBUG
cuda.dependency_type = TYPE_C
QMAKE_EXTRA_COMPILERS += cuda
