#!/bin/sh

export LD_LIBRARY_PATH=/opt/share/gcc/lib64:/opt/share/casa/current/release/lib64
exec /opt/share/casa/current/release/lib64/casapy/bin/python test.py
