
SRCROOT=..
PREFIX=/opt/sakura

SUBDIRS=lib

.PHONY: all lib docs install clean

all: docs

lib:
	(cd lib; make all)

docs: lib
	-mkdir -p $(SRCROOT)/dist/docs/capi
	doxygen Doxyfile

clean:
	(cd lib; make clean)
