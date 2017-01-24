#!/usr/bin/env python3

# To run this program, cd to libsakura/build directory and type as below
# PYTHONPATH=$PWD/python-binding ../test/bench-interpolation-sakura.py 

from bench_template import bench as bench
import libsakurapy
import math
import gc
import numpy
import time

class WSType:
	TYPE = libsakurapy.INTERPOLATION_METHOD_LINEAR # linear interpolation
	ORDER = 1 # effective only when polynomial interpolation
	NUM_BASE = 100
	NUM_ARRAY = 1
	NUM_INTERPOLATED = 1000000 # from bench-interpolation.cc
	INTERP_START = 0.0
	INTERP_END = 101.0
	
class ThreadLocalType:
	def __init__(self, thread_id):
		self.thread_id = thread_id

class InterpolationXBench:
	@staticmethod
	def createWS(id):
		ws = WSType()
		# position list
		base_position = numpy.arange(1, 101, dtype=numpy.float64)
		#print(base_position)
		assert len(base_position) == ws.NUM_BASE
		ws.base_position = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_DOUBLE, base_position)
		# input data
		amp = 10.0
		in_data = base_position + numpy.asarray([amp * ((-1)**i) for i in range(len(base_position))])
		#print(in_data)
		ws.in_data = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_FLOAT, in_data)
		# input mask
		ws.is_valid = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_BOOL, (True,) * ws.NUM_BASE)
		# position list for interpolation
		delta = (ws.INTERP_END - ws.INTERP_START) / float(ws.NUM_INTERPOLATED - 1)
		interp_position = numpy.arange(ws.INTERP_START, ws.INTERP_END + delta / 2, delta, dtype=numpy.float64)
		assert len(interp_position) == ws.NUM_INTERPOLATED
		#print(interp_position)
		ws.interp_position = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_DOUBLE, interp_position)
		# output data
		ws.out_data = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, (ws.NUM_INTERPOLATED, ws.NUM_ARRAY,))
		# output mask
		ws.is_valid_out = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_BOOL, (ws.NUM_INTERPOLATED, ws.NUM_ARRAY,))

		return ws

	@staticmethod
	def destroyWS(id, ws):
		del ws.base_position
		del ws.in_data
		del ws.is_valid
		del ws.interp_position
		del ws.out_data
		del ws.is_valid_out
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
		residual = libsakurapy.interpolate_float_xaxis(
			ws.TYPE, ws.ORDER, 
			ws.NUM_BASE, ws.base_position, ws.NUM_ARRAY, ws.in_data, ws.is_valid,
			ws.NUM_INTERPOLATED, ws.interp_position, ws.out_data, ws.is_valid_out)
		

def main():
	libsakurapy.initialize()
	try:
		bench(InterpolationXBench, 8, 1, 20000)
	finally:
		gc.collect(2)
		libsakurapy.clean_up()

main()
