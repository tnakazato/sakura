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

buf = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, (10, 20, 30, 40, 50))
print(libsakurapy.get_elements_of_aligned_buffer(buf))
del buf

n = 1024*1024*16;
ui8 = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_INT8, (n,))
bl = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_BOOL, (n,))
libsakurapy.uint8_to_bool(n, ui8, bl)
del ui8
del bl

gc.collect(2)

libsakurapy.clean_up()
