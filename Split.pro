TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    findZero.cpp \
    Circuit_F2_1.cpp \
    main.cpp

HEADERS += \
    Circuit_F2_1.h

LIBS += -L/export/home/weixin/lib -laffa

INCLUDEPATH += /export/home/weixin/Desktop/libaffa-0.9.6/src
