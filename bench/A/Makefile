SRCROOT=$(abspath ../..)

CASA=/opt/casapy-34.0.19988-002-64b
SQLITE=/opt/sqlite
LIBPATH=$(SRCROOT)/dist/lib:$(SQLITE)/lib:$(CASA)/lib
CC=gcc
CXX=g++
#CC=icc
#CXX=icc
#LINKER=icpc
LINKER=LD_LIBRARY_PATH=$(LIBPATH) LD_RUN_PATH=$(LIBPATH) $(CXX)
COMMON_FLAGS=-g -DPREFIX='"$(SRCROOT)/dist"' -I$(SRCROOT)/dist/include -I$(SQLITE)/include -I$(CASA)/include -I$(CASA)/include/casacore -Wall -O2 -fopenmp -m64 -march=native #-ftree-vectorize -funroll-loops -floop-interchange -D_XOPEN_SOURCE=600 #-DSQLDEBUG=1 #-D_FILE_OFFSET_BITS=64
CFLAGS=-std=c99 $(COMMON_FLAGS)
CXXFLAGS=$(COMMON_FLAGS)
OBJS=fetchAllPointing.o fetchAllPointing-sqlite.o insertToPointing.o insertToPointing-sqlite.o

vpath %.h       .:$(SRCROOT)/dist/include

.PHONY: all clean
all: fetchAllPointing fetchAllPointing-sqlite insertToPointing insertToPointing-sqlite
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

.c.s:
	$(CC) $(CFLAGS) -S $<

.cc.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

fetchAllPointing.o: fetchAllPointing.cc
fetchAllPointing-sqlite.o: fetchAllPointing-sqlite.cc SQLite.h
insertToPointing.o: insertToPointing.cc
insertToPointing-sqlite.o: insertToPointing-sqlite.cc

clean:
	-rm -f $(OBJS) fetchAllPointing fetchAllPointing-sqlite insertToPointing insertToPointing-sqlite *~

fetchAllPointing: fetchAllPointing.o
	$(LINKER) $(CXXFLAGS) -o $@ $< -L $(SRCROOT)/dist/lib -L $(SQLITE)/lib -L $(CASA)/lib -lSQLiteCDBC -lsqlite3 -lcasacore

fetchAllPointing-sqlite: fetchAllPointing-sqlite.o
	$(LINKER) $(CXXFLAGS) -o $@ $< -L $(SRCROOT)/dist/lib -L $(SQLITE)/lib -L $(CASA)/lib -lSQLiteCDBC -lsqlite3 -lcasacore

insertToPointing: insertToPointing.o
	$(LINKER) $(CXXFLAGS) -o $@ $< -L $(SRCROOT)/dist/lib -L $(SQLITE)/lib -L $(CASA)/lib -lSQLiteCDBC -lsqlite3 -lcasacore

insertToPointing-sqlite: insertToPointing-sqlite.o
	$(LINKER) $(CXXFLAGS) -o $@ $< -L $(SRCROOT)/dist/lib -L $(SQLITE)/lib -L $(CASA)/lib -lSQLiteCDBC -lsqlite3 -lcasacore

test:
	#SAKURA_ROOT=$(SRCROOT)/dist ./msconverter x141.004Ants.ms sqms
	#LD_LIBRARY_PATH=$(SRCROOT)/dist/lib $(SQLITE)/bin/sqlite3 -header -nullvalue '.null.' -init test.sql sqms.tdb
 
