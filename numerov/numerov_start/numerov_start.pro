QT += core
QT -= gui

TARGET = numerov_start
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -O2 -Wall -Wextra
SOURCES += main.cpp

