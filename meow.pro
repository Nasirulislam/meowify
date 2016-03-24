TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -O2 \
    -Wunused-parameter

INCLUDEPATH += $$PWD/3rdParty/SoundTouch

SOURCES += main.cpp \
    pitchshift.cpp \
    dywapitchtrack.c \
    wavfile.cpp \
    3rdParty/SoundTouch/AAFilter.cpp \
    3rdParty/SoundTouch/BPMDetect.cpp \
    3rdParty/SoundTouch/cpu_detect_x86.cpp \
    3rdParty/SoundTouch/FIFOSampleBuffer.cpp \
    3rdParty/SoundTouch/FIRFilter.cpp \
    3rdParty/SoundTouch/InterpolateCubic.cpp \
    3rdParty/SoundTouch/InterpolateLinear.cpp \
    3rdParty/SoundTouch/InterpolateShannon.cpp \
    3rdParty/SoundTouch/PeakFinder.cpp \
    3rdParty/SoundTouch/RateTransposer.cpp \
    3rdParty/SoundTouch/SoundTouch.cpp \
    3rdParty/SoundTouch/sse_optimized.cpp \
    3rdParty/SoundTouch/TDStretch.cpp \
    soundstretch.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    dywapitchtrack.h \
    wavfile.h \
    3rdParty/SoundTouch/AAFilter.h \
    3rdParty/SoundTouch/cpu_detect.h \
    3rdParty/SoundTouch/FIRFilter.h \
    3rdParty/SoundTouch/InterpolateCubic.h \
    3rdParty/SoundTouch/InterpolateLinear.h \
    3rdParty/SoundTouch/InterpolateShannon.h \
    3rdParty/SoundTouch/PeakFinder.h \
    3rdParty/SoundTouch/RateTransposer.h \
    3rdParty/SoundTouch/TDStretch.h \
    3rdParty/SoundTouch/STTypes.h
