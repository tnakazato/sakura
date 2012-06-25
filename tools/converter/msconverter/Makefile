SRCROOT=../../..

CASA=/home/kohji/casahpc2/linux_64b
SQLITE=/opt/sqlite
LIBPATH=$(SRCROOT)/dist/lib:$(SQLITE)/lib:$(CASA)/lib
CC=gcc
CXX=g++
#CC=icc
#CXX=icc
#LINKER=icpc
LINKER=LD_LIBRARY_PATH=$(LIBPATH) LD_RUN_PATH=$(LIBPATH) $(CXX)
COMMON_FLAGS=-g -I$(SRCROOT)/dist/include -I$(SQLITE)/include -I$(CASA)/include/casacore -Wall -O2 -fopenmp -m64 -march=native #-ftree-vectorize -funroll-loops -floop-interchange -D_XOPEN_SOURCE=600 #-DSQLDEBUG=1 #-D_FILE_OFFSET_BITS=64
CFLAGS=-std=c99 $(COMMON_FLAGS)
CXXFLAGS=$(COMMON_FLAGS)
OBJS=msconverter.o

vpath %.h	.:$(SRCROOT)/dist/include

.PHONY: all clean
all: msconverter
	-mkdir $(SRCROOT)/dist/bin
	cp msconverter $(SRCROOT)/dist/bin

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<

.c.s:
	$(CC) $(CFLAGS) -S $<

.cc.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

msconverter.o: msconverter.cc SQLite.h

clean:
	-rm -f $(OBJS) msconverter *~

msconverter: $(OBJS)
	$(LINKER) $(CXXFLAGS) -o $@ $(OBJS) -L $(SRCROOT)/dist/lib -L $(SQLITE)/lib -L $(CASA)/lib -lSQLiteCDBC -lsqlite3 -lcasa_ms

test:	msconverter
	-ln -s /share4hpc/NightlyProfile/PerformanceTest/Data/x141.004Ants.ms .
	-rm sqms[mt].db
	./msconverter -p $(SRCROOT)/dist x141.004Ants.ms sqms
