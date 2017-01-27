#!/usr/bin/env python3

# To run this program, cd to libsakura/build directory and type as below
# PYTHONPATH=$PWD/python-binding ../test/bench-interpolation-sakura.py 

from bench_template import bench as bench
import libsakurapy
import math
import gc
import numpy
import numpy.ma as ma
import time
from scipy import interpolate

class WSType:
	TYPE = libsakurapy.INTERPOLATION_METHOD_LINEAR # linear interpolation
	ORDER = 1 # effective only when polynomial interpolation
	NUM_BASE = 100
	NUM_ARRAY = 1
	NUM_INTERPOLATED = 1000000 # from bench-interpolation.cc, Release
	#NUM_INTERPOLATED = 100000 # Debug
	INTERP_START = 0.0
	INTERP_END = 101.0
	
class ThreadLocalType:
	def __init__(self, thread_id):
		self.thread_id = thread_id

class InterpolationXBench:
	@staticmethod
	def createWS(id):
		working_set = WSType()
		# Known data
		# xi = i - 1
		# yi = xi + amp*(-1)^i
		known_positions = numpy.arange(1, 101, dtype=numpy.float64)
		assert known_positions.size == working_set.NUM_BASE

		amp = 10.0
		offset = numpy.ones(known_positions.size)
		offset[1::2] = -1
		known_values = known_positions + offset*amp

		# Invalid known data marker: no invalid value
		is_invalid = numpy.zeros(known_values.shape,dtype=numpy.bool)

		# Function interpolating known data
		# Note: with current sakuradev01 python3-scipy version 0.13.3, cannot extrapolate like sakura
		# fill_value = fill_like_sakura = (fill_below,fill_above) = (known_values[0],known_values[-1]) : requires scipy >= 0.17.0
		fill_fixed_value = known_values[0]
		marked_known_values = ma.array(known_values,mask=is_invalid) # Mask (=hide) data marked as invalid 
		working_set.interpolating_function = interpolate.interp1d(known_positions,marked_known_values,bounds_error=False,fill_value=fill_fixed_value)

		# Data to interpolate
		delta = (working_set.INTERP_END - working_set.INTERP_START) / float(working_set.NUM_INTERPOLATED - 1)
		wanted_positions = numpy.arange(working_set.INTERP_START, working_set.INTERP_END + delta / 2, delta, dtype=numpy.float64)
		assert wanted_positions.size == working_set.NUM_INTERPOLATED
		working_set.wanted_positions = wanted_positions

		return working_set

	@staticmethod
	def destroyWS(id, ws):
		del ws.wanted_positions
		del ws.interpolating_function
		del ws

	@staticmethod
	def createThreadLocal(thread_id):
		tl = ThreadLocalType(thread_id)
		return tl

	@staticmethod
	def destroyThreadLocal(tread_id, threadLocal):
		del threadLocal

	@staticmethod
	def run(tl, ws, id):
		interpolated_values = ws.interpolating_function(ws.wanted_positions)

def main():
	libsakurapy.initialize()
	try:
		#working_sets, threads, iterations = 8,1,20000 # Release
		working_sets, threads, iterations = 8,8,100 # Debug
		bench(InterpolationXBench, working_sets, threads, iterations)
	finally:
		gc.collect(2)
		libsakurapy.clean_up()

main()
