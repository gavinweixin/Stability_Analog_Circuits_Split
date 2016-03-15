TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    recursiveTaylor.cpp \
    findZero.cpp \
    Circuit_F2_1.cpp \
    main.cpp

HEADERS += \
    Circuit_F2_1.h

LIBS += -L/export/home/weixin/lib -lginac -lcln

INCLUDEPATH += /export/home/weixin/include
