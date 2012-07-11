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
OBJS=get-all-casa.o get-all-sqlite.o

vpath %.h       .:$(SRCROOT)/dist/include

.PHONY: all clean
all: get-all-casa get-all-sqlite

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

.c.s:
	$(CC) $(CFLAGS) -S $<

.cc.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

get-all-casa.o: get-all-casa.cc
get-all-sqlite.o: get-all-sqlite.cc SQLite.h

clean:
	-rm -f $(OBJS) get-all-casa get-all-sqlite *~

get-all-casa: get-all-casa.o
	$(LINKER) $(CXXFLAGS) -o $@ $< -L $(SRCROOT)/dist/lib -L $(SQLITE)/lib -L $(CASA)/lib -lSQLiteCDBC -lsqlite3 -lcasacore

get-all-sqlite: get-all-sqlite.o
	$(LINKER) $(CXXFLAGS) -o $@ $< -L $(SRCROOT)/dist/lib -L $(SQLITE)/lib -L $(CASA)/lib -lSQLiteCDBC -lsqlite3 -lcasacore

test:
	#SAKURA_ROOT=$(SRCROOT)/dist ./msconverter x141.004Ants.ms sqms
	#LD_LIBRARY_PATH=$(SRCROOT)/dist/lib $(SQLITE)/bin/sqlite3 -header -nullvalue '.null.' -init test.sql sqms.tdb
