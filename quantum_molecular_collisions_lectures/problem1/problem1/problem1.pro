QT += core
QT -= gui

TARGET = problem1
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

QMAKE_CXXFLAGS += -std=c++11 -O2 -Wall -Wextra -g

INCLUDEPATH += /usr/local/include/eigen3

SOURCES += main.cpp

