#!/bin/env python

import sys
import time
import gc
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

ctx1D = libsakurapy.create_convolve1D_context(10, libsakurapy.CONVOLVE1D_KERNEL_TYPE_GAUSSIAN, 4, True);
del ctx1D

mask = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_BOOL, (1, 1, 1, 1))
data = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_FLOAT, (1, 2, 3, 4))

result = libsakurapy.compute_statistics(4, data, mask)
print result

del mask
del data

gc.collect(2)

libsakurapy.clean_up()
