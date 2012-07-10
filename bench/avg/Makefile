SQLITE=/opt/sqlite
SAKURA=$(HOME)/sakura/dist
CC=gcc
CXX=g++
#CC=icc
#CXX=icc
#LINKER=icpc
LINKER=LD_LIBRARY_PATH=$(SQLITE)/lib:$(SAKURA)/lib LD_RUN_PATH=$(SQLITE)/lib:$(SAKURA)/lib $(CXX)
COMMON_FLAGS=-g -I$(SQLITE)/include -I$(SAKURA)/include -Wall -O2 -ftree-vectorize -funroll-loops -m64 -march=native -floop-interchange -fopenmp -D_XOPEN_SOURCE=600 #-DSQLDEBUG=1 #-D_FILE_OFFSET_BITS=64
CFLAGS=-std=c99 $(COMMON_FLAGS)
CXXFLAGS=$(COMMON_FLAGS)
OBJS=average.o
SONAME=libaverage-ext.so.1
SOVER=.0.0

.PHONY: all run clean test
all: average $(SONAME)$(SOVER)

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

.c.s:
	$(CC) $(CFLAGS) -S $<

.cc.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(SONAME)$(SOVER): average-ext.cc
	$(CXX) -shared -fPIC $(CXXFLAGS) average-ext.cc -o $(SONAME)$(SOVER) -Wl,-soname,$(SONAME)

average.o: average.cc

clean:
	-rm -f $(OBJS) average *~

average:  average.o
	$(LINKER) $(CXXFLAGS) -o $@ average.o -L$(SAKURA)/lib -L$(SQLITE)/lib -lSQLiteCDBC -lsqlite3

test: all
	-rm out.tdb
	./average x141.001Ants.mdb x141.001Ants.tdb out.tdb
