#!/bin/env python

import sys
import time
from threading import *

sys.path.append(r'.')
import libsakurapy

print dir(libsakurapy);

libsakurapy.initialize()

print libsakurapy.get_current_time()

try:
	libsakurapy.compute_statistics()
except Exception as e:
	print e

mask = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_BOOL, (1, 1, 1, 1))
data = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_FLOAT, (1, 2, 3, 4))

result = libsakurapy.compute_statistics(4, data, mask)
print result

libsakurapy.clean_up()
