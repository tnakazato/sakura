
SRCROOT=../..

SQLITE=/opt/sqlite
CC=gcc
CXX=g++
#CC=icc
#CXX=icc
#LINKER=icpc
LINKER=LD_LIBRARY_PATH=$(SQLITE)/lib LD_RUN_PATH=$(SQLITE)/lib $(CXX)
COMMON_FLAGS=-g -I../include -I$(SQLITE)/include -Wall -O2 -fopenmp -m64 -march=native # -ftree-vectorize -funroll-loops -floop-interchange -D_XOPEN_SOURCE=600 #-DSQLDEBUG=1 #-D_FILE_OFFSET_BITS=64
CFLAGS=-std=c99 $(COMMON_FLAGS)
CXXFLAGS=$(COMMON_FLAGS)
OBJS=SQLite.o
SONAMEwoV=libSQLiteCDBC.so
SONAME=$(SONAMEwoV).1
SOVER=.0.0

vpath %.h ../include

.PHONY: all install clean
all: $(SONAME)$(SOVER)
	-mkdir -p $(SRCROOT)/dist/lib
	cp $(SONAME)$(SOVER) $(SRCROOT)/dist/lib
	-ln -s $(SRCROOT)/dist/lib/$(SONAME)$(SOVER) $(SRCROOT)/dist/lib/$(SONAME)
	-ln -s $(SRCROOT)/dist/lib/$(SONAME) $(SRCROOT)/dist/lib/$(SONAMEwoV)
	-mkdir -p $(SRCROOT)/dist/include
	cp ../include/*.h $(SRCROOT)/dist/include

$(SONAME)$(SOVER): $(OBJS)
	$(CXX) -shared $(CXXFLAGS) $(OBJS) -o $(SONAME)$(SOVER) -Wl,-soname,$(SONAME)

.c.o:
	$(CC) -fPIC $(CFLAGS) -o $@ -c $<

.c.s:
	$(CC) -fPIC $(CFLAGS) -S $<

.cc.o:
	$(CXX) -fPIC $(CXXFLAGS) -o $@ -c $<

SQLite.o: SQLite.cc SQLite.h

clean:
	-rm -f $(OBJS) $(SONAME)$(SOVER) *~
