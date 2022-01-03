#############################################################################
# Makefile for building: StructDynGlass
# Generated by qmake (2.01a) (Qt 4.8.7) on: lun. janv. 3 17:32:28 2022
# Project:  StructDynGlass.pro
# Template: app
# Command: /usr/lib/x86_64-linux-gnu/qt4/bin/qmake -o Makefile StructDynGlass.pro
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
CFLAGS        = -m64 -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -m64 -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
INCPATH       = -I/usr/share/qt4/mkspecs/linux-g++-64 -I. -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I/usr/include/qt4 -I. -I.
LINK          = g++
LFLAGS        = -m64 -Wl,-O1
LIBS          = $(SUBLIBS)  -L/usr/lib/x86_64-linux-gnu -lQtGui -lQtCore -lpthread 
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/lib/x86_64-linux-gnu/qt4/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = defs.cpp \
		dyn_bb.cpp \
		dyn_exp.cpp \
		dyn_isf.cpp \
		eval_isoconf.cpp \
		global.cpp \
		main.cpp \
		nrutil.cpp \
		pbc.cpp \
		read_write_lammps.cpp \
		struct_base.cpp 
OBJECTS       = defs.o \
		dyn_bb.o \
		dyn_exp.o \
		dyn_isf.o \
		eval_isoconf.o \
		global.o \
		main.o \
		nrutil.o \
		pbc.o \
		read_write_lammps.o \
		struct_base.o
DIST          = /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/shared.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		StructDynGlass.pro
QMAKE_TARGET  = StructDynGlass
DESTDIR       = 
TARGET        = StructDynGlass

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: StructDynGlass.pro  /usr/share/qt4/mkspecs/linux-g++-64/qmake.conf /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/shared.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		/usr/lib/x86_64-linux-gnu/libQtGui.prl \
		/usr/lib/x86_64-linux-gnu/libQtCore.prl
	$(QMAKE) -o Makefile StructDynGlass.pro
/usr/share/qt4/mkspecs/common/unix.conf:
/usr/share/qt4/mkspecs/common/linux.conf:
/usr/share/qt4/mkspecs/common/gcc-base.conf:
/usr/share/qt4/mkspecs/common/gcc-base-unix.conf:
/usr/share/qt4/mkspecs/common/g++-base.conf:
/usr/share/qt4/mkspecs/common/g++-unix.conf:
/usr/share/qt4/mkspecs/qconfig.pri:
/usr/share/qt4/mkspecs/features/qt_functions.prf:
/usr/share/qt4/mkspecs/features/qt_config.prf:
/usr/share/qt4/mkspecs/features/exclusive_builds.prf:
/usr/share/qt4/mkspecs/features/default_pre.prf:
/usr/share/qt4/mkspecs/features/release.prf:
/usr/share/qt4/mkspecs/features/default_post.prf:
/usr/share/qt4/mkspecs/features/shared.prf:
/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf:
/usr/share/qt4/mkspecs/features/warn_on.prf:
/usr/share/qt4/mkspecs/features/qt.prf:
/usr/share/qt4/mkspecs/features/unix/thread.prf:
/usr/share/qt4/mkspecs/features/moc.prf:
/usr/share/qt4/mkspecs/features/resources.prf:
/usr/share/qt4/mkspecs/features/uic.prf:
/usr/share/qt4/mkspecs/features/yacc.prf:
/usr/share/qt4/mkspecs/features/lex.prf:
/usr/share/qt4/mkspecs/features/include_source_dir.prf:
/usr/lib/x86_64-linux-gnu/libQtGui.prl:
/usr/lib/x86_64-linux-gnu/libQtCore.prl:
qmake:  FORCE
	@$(QMAKE) -o Makefile StructDynGlass.pro

dist: 
	@$(CHK_DIR_EXISTS) .tmp/StructDynGlass1.0.0 || $(MKDIR) .tmp/StructDynGlass1.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/StructDynGlass1.0.0/ && $(COPY_FILE) --parents defs.h dyn_bb.h dyn_exp.h dyn_isf.h eval_isoconf.h global.h nrutil.h pbc.h read_write_lammps.h struct_base.h .tmp/StructDynGlass1.0.0/ && $(COPY_FILE) --parents defs.cpp dyn_bb.cpp dyn_exp.cpp dyn_isf.cpp eval_isoconf.cpp global.cpp main.cpp nrutil.cpp pbc.cpp read_write_lammps.cpp struct_base.cpp .tmp/StructDynGlass1.0.0/ && (cd `dirname .tmp/StructDynGlass1.0.0` && $(TAR) StructDynGlass1.0.0.tar StructDynGlass1.0.0 && $(COMPRESS) StructDynGlass1.0.0.tar) && $(MOVE) `dirname .tmp/StructDynGlass1.0.0`/StructDynGlass1.0.0.tar.gz . && $(DEL_FILE) -r .tmp/StructDynGlass1.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


check: first

mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_moc_header_make_all:
compiler_moc_header_clean:
compiler_rcc_make_all:
compiler_rcc_clean:
compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_uic_make_all:
compiler_uic_clean:
compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile

defs.o: defs.cpp defs.h \
		nrutil.h \
		read_write_lammps.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o defs.o defs.cpp

dyn_bb.o: dyn_bb.cpp dyn_bb.h \
		defs.h \
		nrutil.h \
		read_write_lammps.h \
		pbc.h \
		eval_isoconf.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dyn_bb.o dyn_bb.cpp

dyn_exp.o: dyn_exp.cpp dyn_exp.h \
		defs.h \
		nrutil.h \
		read_write_lammps.h \
		pbc.h \
		eval_isoconf.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dyn_exp.o dyn_exp.cpp

dyn_isf.o: dyn_isf.cpp dyn_isf.h \
		defs.h \
		nrutil.h \
		read_write_lammps.h \
		pbc.h \
		eval_isoconf.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o dyn_isf.o dyn_isf.cpp

eval_isoconf.o: eval_isoconf.cpp eval_isoconf.h \
		defs.h \
		nrutil.h \
		read_write_lammps.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o eval_isoconf.o eval_isoconf.cpp

global.o: global.cpp global.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o global.o global.cpp

main.o: main.cpp defs.h \
		nrutil.h \
		read_write_lammps.h \
		pbc.h \
		dyn_bb.h \
		dyn_exp.h \
		dyn_isf.h \
		struct_base.h \
		global.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

nrutil.o: nrutil.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o nrutil.o nrutil.cpp

pbc.o: pbc.cpp pbc.h \
		defs.h \
		nrutil.h \
		read_write_lammps.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o pbc.o pbc.cpp

read_write_lammps.o: read_write_lammps.cpp defs.h \
		nrutil.h \
		read_write_lammps.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o read_write_lammps.o read_write_lammps.cpp

struct_base.o: struct_base.cpp struct_base.h \
		defs.h \
		nrutil.h \
		read_write_lammps.h \
		dyn_bb.h \
		pbc.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o struct_base.o struct_base.cpp

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

