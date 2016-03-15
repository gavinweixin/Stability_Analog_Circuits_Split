TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    findZero.cpp \
    Circuit_F2_1.cpp \
    main.cpp \
    alterOutput.cpp \
    Circuit.cpp

HEADERS += \
    Circuit_F2_1.h \
    Circuit.h

LIBS += -L/export/home/weixin/lib -laffa

INCLUDEPATH += /export/home/weixin/include
