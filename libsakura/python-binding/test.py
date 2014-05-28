#!/bin/env python

import sys
import time
import math
import gc
from threading import *

import libsakurapy

print dir(libsakurapy);

libsakurapy.initialize()

print libsakurapy.get_current_time()

def test_stats():
	try:
		libsakurapy.compute_statistics()
	except Exception as e:
		print e

	mask = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_BOOL, (1, 1, 1, 1))
	data = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_FLOAT, (1, 2, 3, 4))
	result = libsakurapy.compute_statistics(4, data, mask)
	print result
	del mask
	del data

def test_grid():
	def total_elements(dim):
		return reduce(lambda x, y: x * y, dim, 1)

	num_spectra = 1024*64
	start_spectrum = 0
	end_spectrum = num_spectra
	spectrum_mask = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_BOOL, (num_spectra,))
	x = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_DOUBLE, (num_spectra,))
	y = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_DOUBLE, (num_spectra,))
	support = 10
	sampling = 2
	num_polarizations = 4
	num_polarizations_for_grid = 2
	num_channels = 2048
	num_channels_for_grid = 1024
	polarization_map = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_INT32,
		map(lambda x: x % num_polarizations_for_grid, range(num_polarizations)))
	channel_map = libsakurapy.new_aligned_buffer(libsakurapy.TYPE_INT32,
		map(lambda x: x % num_channels_for_grid, range(num_channels)))

	dim = (num_spectra, num_polarizations, num_channels)
	print total_elements(dim)
	mask = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_INT8, dim)
	mask = libsakurapy.uint8_to_bool(total_elements(dim),
		mask, libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_BOOL, dim))
	value = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, dim)
	weight = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, (num_spectra, num_channels))
	num_convolution_table = int(math.ceil(math.sqrt(2.)*(support+1)*sampling))
	convolution_table = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, (num_convolution_table, ))

	width, height = (160, 100)
	weight_sum = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_DOUBLE, (num_polarizations_for_grid, num_channels_for_grid))
	dim = (height, width, num_polarizations_for_grid, num_channels_for_grid)
	weight_of_grid = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, dim)
	grid = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, dim)

	libsakurapy.grid_convolving(num_spectra, start_spectrum, end_spectrum,
		spectrum_mask, x, y,
		support, sampling,
		num_polarizations, polarization_map,
		num_channels, channel_map,
		mask, value, weight, False,
		num_convolution_table, convolution_table,
		num_polarizations_for_grid, num_channels_for_grid,
		width, height,
		weight_sum, weight_of_grid, grid)

def test_AB():
	buf = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, (10, 20, 30, 40, 50))
	print(libsakurapy.get_elements_of_aligned_buffer(buf))
	del buf

def test_logical():
	n = 1024*1024*16
	dim = (n,)
	
	buf = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, dim)
	bl = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_BOOL, dim)
	libsakurapy.set_false_float_if_nan_or_inf(n, buf, bl)
	del buf
	del bl

	ui8 = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_INT8, dim)
	bl = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_BOOL, dim)
	libsakurapy.uint8_to_bool(n, ui8, bl)
	ui32 = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_INT32, dim)
	libsakurapy.uint32_to_bool(n, ui32, bl)
	bl2 = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_BOOL, dim)
	libsakurapy.invert_bool(n, bl, bl2)
	libsakurapy.invert_bool(n, bl, bl) # in place
	del ui8
	del ui32
	del bl
	del bl2

def test_convolve1D():
	ctx1D = libsakurapy.create_convolve1D_context(10, libsakurapy.CONVOLVE1D_KERNEL_TYPE_GAUSSIAN, 4, True);
	del ctx1D

def testAll():
	test_AB()
	test_stats()
	test_grid()
	test_logical()
	test_convolve1D()

testAll()
gc.collect(2)

libsakurapy.clean_up()
