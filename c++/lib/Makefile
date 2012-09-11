
SRCROOT=../..

SQLITE=/opt/sqlite
CC=gcc
CXX=g++
#CC=icc
#CXX=icc
#LINKER=icpc
LINKER=LD_LIBRARY_PATH=$(SQLITE)/lib LD_RUN_PATH=$(SQLITE)/lib $(CXX)
COMMON_FLAGS=-g -O2 -m64 -march=native -Wall -fopenmp
COMMON_FLAGS+=-I../include -I$(SQLITE)/include
COMMON_FLAGS+=-D_XOPEN_SOURCE=600 # -ftree-vectorize -funroll-loops -floop-interchange #-DSQLDEBUG=1 #-D_FILE_OFFSET_BITS=64
CFLAGS=-std=c99 $(COMMON_FLAGS)
CXXFLAGS=$(COMMON_FLAGS)
OBJS=SQLite.o
SONAME_CDBCwoV=libSQLiteCDBC.so
SONAME_CDBC=$(SONAME_CDBCwoV).1
SOVER_CDBC=.0.0

SONAME_SQLFuncwoV=libSQLiteFunc.so
SONAME_SQLFunc=$(SONAME_SQLFuncwoV).1
SOVER_SQLFunc=.0.0

SONAME_MMapModwoV=libSQLiteMMapVTable.so
SONAME_MMapMod=$(SONAME_MMapModwoV).1
SOVER_MMapMod=.0.0

vpath %.h ../include

.PHONY: all install clean
all: $(SONAME_CDBC)$(SOVER_CDBC) $(SONAME_SQLFunc)$(SOVER_SQLFunc) $(SONAME_MMapMod)$(SOVER_MMapMod)
	-mkdir -p $(SRCROOT)/dist/lib
	cp $(SONAME_CDBC)$(SOVER_CDBC) $(SRCROOT)/dist/lib
	-ln -s $(SRCROOT)/dist/lib/$(SONAME_CDBC)$(SOVER_CDBC) $(SRCROOT)/dist/lib/$(SONAME_CDBC)
	-ln -s $(SRCROOT)/dist/lib/$(SONAME_CDBC) $(SRCROOT)/dist/lib/$(SONAME_CDBCwoV)
	cp $(SONAME_SQLFunc)$(SOVER_SQLFunc) $(SRCROOT)/dist/lib
	-ln -s $(SRCROOT)/dist/lib/$(SONAME_SQLFunc)$(SOVER_SQLFunc) $(SRCROOT)/dist/lib/$(SONAME_SQLFunc)
	-ln -s $(SRCROOT)/dist/lib/$(SONAME_SQLFunc) $(SRCROOT)/dist/lib/$(SONAME_SQLFuncwoV)
	cp $(SONAME_MMapMod)$(SOVER_MMapMod) $(SRCROOT)/dist/lib
	-ln -s $(SRCROOT)/dist/lib/$(SONAME_MMapMod)$(SOVER_MMapMod) $(SRCROOT)/dist/lib/$(SONAME_MMapMod)
	-ln -s $(SRCROOT)/dist/lib/$(SONAME_MMapMod) $(SRCROOT)/dist/lib/$(SONAME_MMapModwoV)
	-mkdir -p $(SRCROOT)/dist/include
	cp ../include/*.h $(SRCROOT)/dist/include

$(SONAME_CDBC)$(SOVER_CDBC): SQLite.o
	$(CXX) -shared $(CXXFLAGS) SQLite.o -o $(SONAME_CDBC)$(SOVER_CDBC) -Wl,-soname,$(SONAME_CDBC)

$(SONAME_SQLFunc)$(SOVER_SQLFunc): sql-ext.o
	$(CXX) -shared $(CXXFLAGS) sql-ext.o -o $(SONAME_SQLFunc)$(SOVER_SQLFunc) -Wl,-soname,$(SONAME_SQLFunc)

$(SONAME_MMapMod)$(SOVER_MMapMod): libSQLiteMMapVTable.o
	$(CXX) -shared $(CXXFLAGS) libSQLiteMMapVTable.o -o $(SONAME_MMapMod)$(SOVER_MMapMod) -Wl,-soname,$(SONAME_MMapMod)

.c.o:
	$(CC) -fPIC $(CFLAGS) -o $@ -c $<

.c.s:
	$(CC) -fPIC $(CFLAGS) -S $<

.cc.o:
	$(CXX) -fPIC $(CXXFLAGS) -o $@ -c $<

SQLite.o: SQLite.cc SQLite.h

sql-ext.o: sql-ext.cc

libSQLiteMMapVTable.o: libSQLiteMMapVTable.cc

clean:
	-rm -f $(OBJS) $(SONAME_CDBC)$(SOVER_CDBC) $(SONAME_SQLFunc)$(SOVER_SQLFunc) $(SONAME_MMapMod)$(SOVER_MMapMod) *~
