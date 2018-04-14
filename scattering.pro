QT += core
QT -= gui

TARGET = scattering
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -lm -O2 -Wall -Wextra

SOURCES += main.cpp

