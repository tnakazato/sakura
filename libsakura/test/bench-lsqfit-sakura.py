#!/usr/bin/env python3

# To run this program, cd to libsakura/build directory and type as below
# PYTHONPATH=$PWD/python-binding ../test/bench-lsqfit-sakura.py 

from bench_template import bench as bench
import libsakurapy
import math
import gc
import numpy
import time

class WSType:
	ORDER = 14
	NUM_COEFF = ORDER + 1
	N = 1000000

class ThreadLocalType:
	def __init__(self, thread_id):
		self.thread_id = thread_id

class MyBench:
	@staticmethod
	def createWS(id):
		ws = WSType()
		ws.is_valid = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_BOOL, (True,) * ws.N)
		ws.coeff = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_DOUBLE, (ws.NUM_COEFF,))
		ws.is_valid_out = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_BOOL, (ws.N,))

		period1 = ws.N
		amplitude1 = 5.0

		period2 = 0.695 * ws.N
		amplitude2 = 5.0
		x = numpy.arange(ws.N, dtype=float)
		y = amplitude1 * numpy.sin(2.0 * numpy.pi * x / period1) \
		    + amplitude2 * numpy.cos(2.0 * numpy.pi * x / period2)
		edge = 0.1 * ws.N
		y0 = y[:edge].mean()
		y1 = y[-edge:].mean()
		x0 = x[0]
		x1 = x[-1]
		l = y0 + (y1 - y0) / (x1 - x0) * x
		z = y - l
		ws.data = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_FLOAT, z)
		return ws

	@staticmethod
	def destroyWS(id, ws):
		del ws.data
		del ws.is_valid
		del ws.coeff
		del ws.is_valid_out
		del ws

	@staticmethod
	def createThreadLocal(thread_id):
		tl = ThreadLocalType(thread_id)
		tl.context = libsakurapy.create_baseline_context(libsakurapy.BASELINE_TYPE_POLYNOMIAL, WSType.ORDER, WSType.N)
		return tl

	@staticmethod
	def destroyThreadLocal(tread_id, threadLocal):
		del threadLocal

	@staticmethod
	def run(tl, ws, id):
		residual = libsakurapy.lsqfit_polynomial(tl.context, WSType.ORDER,
				WSType.N, ws.data, ws.is_valid, 0.0001, 1,
				WSType.NUM_COEFF, ws.coeff, None, None,
				ws.is_valid_out)


def main():
	libsakurapy.initialize()
	try:
		bench(MyBench, 8, 1, 20000)
	finally:
		gc.collect(2)
		libsakurapy.clean_up()

main()
