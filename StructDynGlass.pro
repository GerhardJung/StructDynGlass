######################################################################
# Automatically generated by qmake (3.1) Thu Apr 21 08:47:23 2022
######################################################################

TEMPLATE = app
TARGET = StructDynGlass
INCLUDEPATH += .

# The following define makes your compiler warn you if you use any
# feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input
HEADERS += defs.h \
           dyn_bb.h \
           dyn_exp.h \
           dyn_isf.h \
           dyn_msd.h \
           dyn_rearrange_patinet.h \
           eval_isoconf.h \
           eval_struct.h \
           global.h \
           nrutil.h \
           pbc.h \
           read_write_lammps.h \
           struct_base.h \
           struct_filion.h \
           struct_gnn.h \
           struct_ml.h \
           struct_read.h \
           struct_soft_modes.h \
           struct_voronoi.h
SOURCES += defs.cpp \
           dyn_bb.cpp \
           dyn_exp.cpp \
           dyn_isf.cpp \
           dyn_msd.cpp \
           dyn_rearrange_patinet.cpp \
           eval_isoconf.cpp \
           eval_struct.cpp \
           global.cpp \
           main.cpp \
           nrutil.cpp \
           pbc.cpp \
           read_write_lammps.cpp \
           struct_base.cpp \
           struct_filion.cpp \
           struct_gnn.cpp \
           struct_ml.cpp \
           struct_read.cpp \
           struct_soft_modes.cpp \
           struct_voronoi.cpp
