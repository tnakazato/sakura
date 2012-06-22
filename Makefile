
SRCROOT=.

SUBDIRS=sql c++/lib tools/converter/msconverter

.PHONY: all install clean

all:
	for dir in $(SUBDIRS) ; do (cd $$dir; make all); done

clean:
	for dir in $(SUBDIRS) ; do (cd $$dir; make clean); done
	-rm -rf dist
