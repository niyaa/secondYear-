TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    startup1d.cpp \
    codes1d.cpp \
    jacobip.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    startup1d.h \
    all.h \
    codes1d.h

