#!/usr/bin/env python3

# To run this program, cd to libsakura/build directory and type as below
# PYTHONPATH=$PWD/python-binding ../test/bench-lsqfit-sakura.py 

from bench_template import bench as bench
import gc
import numpy
import time

class WSType:
	ORDER = 14
	NUM_COEFF = ORDER + 1
	N = 100#0000

class ThreadLocalType:
	def __init__(self, thread_id):
		self.thread_id = thread_id

class LSQFitBench:
	@staticmethod
	def createWS(id):
		ws = WSType()

		period1 = ws.N
		amplitude1 = 5.0

		period2 = 0.695 * ws.N
		amplitude2 = 5.0
		ws.x_data = numpy.arange(ws.N, dtype=float)
		y = amplitude1 * numpy.sin(2.0 * numpy.pi * ws.x_data / period1) \
		    + amplitude2 * numpy.cos(2.0 * numpy.pi * ws.x_data / period2)
		edge = int(0.1 * ws.N)
		y0 = y[:edge].mean()
		y1 = y[-edge:].mean()
		x0 = ws.x_data[0]
		x1 = ws.x_data[-1]
		l = y0 + (y1 - y0) / (x1 - x0) * ws.x_data
		z = y - l
		ws.masked_data = numpy.ma.masked_array(z, mask=False, dtype=float, copy=True)
		return ws

	@staticmethod
	def destroyWS(id, ws):
		del ws.masked_data
		del ws.x_data
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
		coeff = numpy.ma.polyfit(ws.x_data, ws.masked_data, WSType.ORDER, full=False, cov=False)


def main():
	try:
		bench(LSQFitBench, 8, 1, 20000)
	finally:
		gc.collect(2)

main()
