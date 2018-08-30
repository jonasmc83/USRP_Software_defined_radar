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

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = USRP_Pulse_Doppler_Radar
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    plot.cpp \
    processing.cpp \
    simulatedtarget.cpp \
    tracker.cpp \
    ppidialog.cpp \
    detection.cpp \
    cameradialog.cpp \
    waveform.cpp

HEADERS  += mainwindow.h \
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

FORMS    += mainwindow.ui \
    ppidialog.ui \
    cameradialog.ui

LIBS += -lqwt -L/usr/local/lib -luhd -lboost_filesystem -lboost_system -lboost_thread -lboost_chrono -lgomp -lopencv_core -lopencv_highgui -lopencv_video -lopencv_imgproc

INCLUDEPATH = /usr/local/include

DEFINES += QWT_DLL

CONFIG += qwt

QMAKE_CXXFLAGS += -std=c++11 -fopenmp -Wno-sign-compare

OTHER_FILES += \
    processing_gpu.cu \
    processing_gpu_double.cu

# CUDA Specific stuff
#QMAKE_CC = nvcc
#QMAKE_CXX = nvcc

CUDA_SOURCES += processing_gpu.cu processing_gpu_double.cu
CUDA_DIR = /usr/local/cuda-8.0            # Path to cuda toolkit install
CUDA_OBJECTS_DIR = .
SYSTEM_NAME = x86_64         # Depending on your system either 'Win32', 'x64', or 'Win64'
SYSTEM_TYPE = 64            # '32' or '64', depending on your system
CUDA_ARCH = sm_30           # Type of CUDA architecture, for example 'compute_10', 'compute_11', 'sm_10'
NVCC_OPTIONS = --use_fast_math

# include paths
INCLUDEPATH += $$CUDA_DIR/include

# library directories
QMAKE_LIBDIR += $$CUDA_DIR/lib64
# Add the necessary libraries
LIBS += -lcuda -lcudart -lcufft

# The following makes sure all path names (which often include spaces) are put between quotation marks
CUDA_INC = $$join(INCLUDEPATH,'" -I"','-I"','"')

# Configuration of the Cuda compiler
CONFIG(debug, debug|release) {
    # Debug mode
    cuda_d.input = CUDA_SOURCES
    cuda_d.output = $$CUDA_OBJECTS_DIR/${QMAKE_FILE_BASE}_cuda.o
    cuda_d.commands = $$CUDA_DIR/bin/nvcc -D_DEBUG $$NVCC_OPTIONS $$CUDA_INC $$LIBS --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
    cuda_d.dependency_type = TYPE_C
    QMAKE_EXTRA_COMPILERS += cuda_d
}
else {
    # Release mode
    cuda.input = CUDA_SOURCES
    cuda.output = $$CUDA_OBJECTS_DIR/${QMAKE_FILE_BASE}_cuda.o
    cuda.commands = $$CUDA_DIR/bin/nvcc $$NVCC_OPTIONS $$CUDA_INC $$LIBS --machine $$SYSTEM_TYPE -arch=$$CUDA_ARCH -c -o ${QMAKE_FILE_OUT} ${QMAKE_FILE_NAME}
    cuda.dependency_type = TYPE_C
    QMAKE_EXTRA_COMPILERS += cuda
}
