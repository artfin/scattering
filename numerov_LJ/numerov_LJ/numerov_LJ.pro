QT += core
QT -= gui

TARGET = numerov_LJ
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -Wall -Wextra -O2
SOURCES += main.cpp

