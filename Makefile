
SRCROOT=.

SUBDIRS=c++/lib tools/converter/msconverter

.PHONY: all install clean

all:
	for dir in $(SUBDIRS) ; do (cd $$dir; make all); done

clean:
	-rm -rf dist
