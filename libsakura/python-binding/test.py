#!/bin/env python

"""
# @SAKURA_LICENSE_HEADER_START@
# Copyright (C) 2013-2016
# National Astronomical Observatory of Japan
# 2-21-1, Osawa, Mitaka, Tokyo, 181-8588, Japan.
# 
# This file is part of Sakura.
# 
# Sakura is free software: you can redistribute it and/or modify it under 
# the terms of the GNU Lesser General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your 
# option) any later version.
# 
# Sakura is distributed in the hope that it will be useful, but WITHOUT 
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public License 
# along with Sakura.  If not, see <http://www.gnu.org/licenses/>.
# @SAKURA_LICENSE_HEADER_END@
"""

import sys
import time
import math
import gc
import numpy
from threading import *

import libsakurapy

print dir(libsakurapy);

libsakurapy.initialize()

print time.time()

def test_stats():
	try:
		num = 4
		mask = numpy.zeros(num, dtype=numpy.bool)
		mask[:] = True
		data = numpy.arange(4, dtype=numpy.float32)
		libsakurapy.compute_statistics(num, data, mask)
	except Exception as e:
		print e

	num = 4
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (num,))
	mask[:] = True
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (num,))
	data[:] = range(num)
	result = libsakurapy.compute_statistics(4, data, mask)
	print result
	del mask
	del data
	
def test_mad():
	num0 = 7
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (num0,))
	data[:] = range(num0)
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (num0,))
	mask[:5] = True
	mask[5:] = False
	num = libsakurapy.sort_data_densely(num0, mask, data)
	print num0, num
	#num = 5
	data_in = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (num,))
	data_in[:] = data[:5]
	print data_in
	data_out = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (num,))
	result = libsakurapy.compute_mad(num, data_in, data_out)
	print result
	
	del data
	del mask
	del data_in
	del data_out
	del result

def test_grid():
	def total_elements(dim):
		return reduce(lambda x, y: x * y, dim, 1)

	num_spectra = 1024*64
	start_spectrum = 0
	end_spectrum = num_spectra
	spectrum_mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (num_spectra,))
	x = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_DOUBLE, (num_spectra,))
	y = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_DOUBLE, (num_spectra,))
	support = 10
	sampling = 2
	num_polarizations = 4
	num_polarizations_for_grid = 2
	num_channels = 2048
	num_channels_for_grid = 1024
	polarization_map = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, (num_polarizations,))
	polarization_map[:] = map(lambda x: x % num_polarizations_for_grid, range(num_polarizations))
	channel_map = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, (num_channels,))
	channel_map[:] = map(lambda x: x % num_channels_for_grid, range(num_channels))

	dim = (num_spectra, num_polarizations, num_channels)
	print total_elements(dim)
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT8, dim)
	mask = libsakurapy.uint8_to_bool(total_elements(dim),
		mask, libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, dim))
	value = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, dim)
	weight = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (num_spectra, num_channels))
	num_convolution_table = int(math.ceil(math.sqrt(2.)*(support+1)*sampling))
	convolution_table = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (num_convolution_table, ))

	width, height = (160, 100)
	weight_sum = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_DOUBLE, (num_polarizations_for_grid, num_channels_for_grid))
	dim = (height, width, num_polarizations_for_grid, num_channels_for_grid)
	weight_of_grid = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, dim)
	grid = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, dim)

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

# def test_AB():
# 	buf = libsakurapy.new_uninitialized_aligned_buffer(libsakurapy.TYPE_FLOAT, (10, 20, 30, 40, 50))
# 	print(libsakurapy.get_elements_of_aligned_buffer(buf))
# 	del buf

def test_logical():
	n = 1024*1024*16
	dim = (n,)
	
	buf = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, dim)
	buf[0] = 1.
	buf[1] = numpy.nan
	bl = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, dim)
	libsakurapy.set_false_float_if_nan_or_inf(n, buf, bl)
	print 'set_false_float_if_nan_or_inf: data {0} mask {1}'.format(buf[:2], bl[:2])
	del buf
	del bl

	ui8 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT8, dim)
	ui8[0] = 0
	ui8[1] = 1
	bl = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, dim)
	libsakurapy.uint8_to_bool(n, ui8, bl)
	print 'uint8_to_bool: before {0} after {1}'.format(ui8[:2], bl[:2])

	ui32 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, dim)
	ui32[0] = 0
	ui32[1] = 1
	libsakurapy.uint32_to_bool(n, ui32, bl)
	print 'uint32_to_bool: before {0} after {1}'.format(ui32[:2], bl[:2])
	
	arr = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, dim)
	arr[:] = True
	arr2 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, dim)
	libsakurapy.invert_bool(n, arr, arr2)
	print 'invert_bool: before {0} after {1}'.format(arr[0], arr2[0])
	arr_save = arr[0]
	libsakurapy.invert_bool(n, arr, arr) # in place
	print 'invert_bool (in-place): before {0} after {1}'.format(arr_save, arr[0])

	del ui8
	del ui32
	del bl
	del arr
	del arr2

	src1 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (4,))
	src1[::2] = True
	src1[1::2] = False
	src2 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (4,))
	src2[:2] = True
	src2[2:] = False
	dst = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (4,))
	out = libsakurapy.logical_and(4, src1, src2, dst)
	assert out is dst
	print 'logical_and: \n\tsrc1 {0}\n\tsrc2 {1}\n\tdst {2}'.format(src1, src2, out)
	
def test_range():
	# Test set_true_int_in_ranges_exclusive
	n = 1024*1024*16
	dim = (n,)
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, dim)
	data[:5] = [0, 1, 2, 205, 300]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, dim)
	lower = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, (3,))
	lower[:] = (1, 100, 200)
	upper = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, (3,))
	upper[:] = (10, 110, 210)
	result = libsakurapy.set_true_int_in_ranges_exclusive(n, data, 3, lower, upper, mask)
	print 'set_true_int_in_ranges_exclusive: data {0} mask {1}'.format(data[:5], mask[:5])
	del n, dim, data, mask, lower, upper, result

	# Test set_true_float_in_ranges_exclusive
	ndata = 4
	dataf = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	dataf[:] = [0., 2., 1., 3.]
	maskf = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	maskf[:] = [True, True, True, True]
	lowerf = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (1,))
	lowerf[:] = [0.5]
	upperf = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (1,))
	upperf[:] = [2.5]
	result = libsakurapy.set_true_float_in_ranges_exclusive(ndata,dataf,1,lowerf,upperf,maskf)
	print 'set_true_float_ranges_exclusive: data {0} mask {1}'.format(dataf, maskf)
	del dataf, ndata, maskf, lowerf, upperf, result
	
	# Test set_true_float_if_greater_than
	ndata = 4
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	data[:] = [0., 2., 1., 3.]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True, True, True, True]
	threshold = numpy.float32(1.0)
	result = libsakurapy.set_true_float_if_greater_than(ndata, data, threshold, mask)
	print 'set_true_float_if_greater_than: threshold {0} data {1} mask {2}'.format(threshold, data, mask)
	del ndata, data, mask, threshold
	
	# Test set_true_float_if_greater_than_or_equal
	ndata = 4
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	data[:] = [0., 2., 1., 3.]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True, True, True, True]
	threshold = numpy.float32(1.0)
	result = libsakurapy.set_true_float_if_greater_than_or_equal(ndata, data, threshold, mask)
	print 'set_true_float_if_greater_than_or_equal: threshold {0} data {1} mask {2}'.format(threshold, data, mask)
	del ndata, data, mask, threshold
	
	# Test set_true_float_if_less_than
	ndata = 4
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	data[:] = [0., 2., 1., 3.]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True, True, True, True]
	threshold = numpy.float32(1.0)
	result = libsakurapy.set_true_float_if_less_than(ndata, data, threshold, mask)
	print 'set_true_float_if_less_than: threshold {0} data {1} mask {2}'.format(threshold, data, mask)
	del ndata, data, mask, threshold
	
	# Test set_true_float_if_less_than_or_equal
	ndata = 4
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	data[:] = [0., 2., 1., 3.]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True, True, True, True]
	threshold = numpy.float32(1.0)
	result = libsakurapy.set_true_float_if_less_than_or_equal(ndata, data, threshold, mask)
	print 'set_true_float_if_less_than_or_equal: threshold {0} data {1} mask {2}'.format(threshold, data, mask)
	del ndata, data, mask, threshold
	
	# Test set_true_int_if_greater_than
	ndata = 4
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, (ndata,))
	data[:] = [0, 2, 1, 3]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True, True, True, True]
	threshold = numpy.int32(1)
	result = libsakurapy.set_true_int_if_greater_than(ndata, data, threshold, mask)
	print 'set_true_int_if_greater_than: threshold {0} data {1} mask {2}'.format(threshold, data, mask)
	del ndata, data, mask, threshold
	
	# Test set_true_int_if_greater_than_or_equal
	ndata = 4
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, (ndata,))
	data[:] = [0, 2, 1, 3]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True, True, True, True]
	threshold = numpy.int32(1)
	result = libsakurapy.set_true_int_if_greater_than_or_equal(ndata, data, threshold, mask)
	print 'set_true_int_if_greater_than_or_equal: threshold {0} data {1} mask {2}'.format(threshold, data, mask)
	del ndata, data, mask, threshold
	
	# Test set_true_int_if_less_than
	ndata = 4
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, (ndata,))
	data[:] = [0, 2, 1, 3]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True, True, True, True]
	threshold = numpy.int32(1)
	result = libsakurapy.set_true_int_if_less_than(ndata, data, threshold, mask)
	print 'set_true_int_if_less_than: threshold {0} data {1} mask {2}'.format(threshold, data, mask)
	del ndata, data, mask, threshold
	
	# Test set_true_int_if_less_than_or_equal
	ndata = 4
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_INT32, (ndata,))
	data[:] = [0, 2, 1, 3]
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True, True, True, True]
	threshold = numpy.int32(1)
	result = libsakurapy.set_true_int_if_less_than_or_equal(ndata, data, threshold, mask)
	print 'set_true_int_if_less_than_or_equal: threshold {0} data {1} mask {2}'.format(threshold, data, mask)
	del ndata, data, mask, threshold
	
def test_bit():
	# Test operate_bits_uint8_not
	ndata = 8	
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True]*4+[False]*4
	data8 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_UINT8, (ndata,))
	data8[:] = [0, 2, 1, 3]*2
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_UINT8, (ndata,))
	out = libsakurapy.operate_bits_uint8_not(ndata, data8, mask, result)
	print 'operate_bits_unit8_not: before {0} after {1}'.format(data8, result)
	del ndata, mask, data8, result, out
	
	# Test operate_bits_uint32_not
	ndata = 8	
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True]*4+[False]*4
	data32 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_UINT32, (ndata,))
	data32[:] = [0, 2, 1, 3]*2
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_UINT32, (ndata,))
	out = libsakurapy.operate_bits_uint32_not(ndata, data32, mask, result)
	print 'operate_bits_unit32_not: before {0} after {1}'.format(data32, result)
	del ndata, mask, data32, result, out

	# Test operate_bits_uint8_or
	ndata = 8	
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True]*4+[False]*4
	data8 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_UINT8, (ndata,))
	data8[:] = [0, 2, 1, 3]*2
	print 'data type for data8: {0}'.format(data8.dtype)
	#ndata = libsakurapy.get_elements_of_aligned_buffer(data8)[0]
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_UINT8, (ndata,))
	out = libsakurapy.operate_bits_uint8_or(2,ndata,data8,mask,result)
	del mask, data8, result, ndata, out

	# Test operate_bits_uint32_or
	ndata = 8
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = [True]*4+[False]*4
	data32 = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_UINT32, (ndata,))
	data32[:] = [0, 2, 1, 3]*2
	#ndata = libsakurapy.get_elements_of_aligned_buffer(data32)[0]
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_UINT32, (ndata,))
	out = libsakurapy.operate_bits_uint32_or(2,ndata,data32,mask,result)
	del mask, data32, result, ndata, out

def test_interpolate():
	# interpolate in Y-axis
	nchan = 4
	yin = [0., 1.]
	yout = [0.75]
	zin = [float(6.0)]*nchan + [float(5.0)]*nchan
	nbase = len(yin)
	npos = len(yout)
	order = 1
	zindata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (nbase, nchan,))
	zindata[0] = zin[:nchan]
	zindata[1] = zin[nchan:]
	yindata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_DOUBLE, (nbase,))
	yindata[:] = yin
	youtdata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_DOUBLE, (npos,))
	youtdata[:] = yout
	zoutdata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (npos, nchan,))
	inmask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (nbase, nchan,))
	inmask[:] = True
	outmask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (npos, nchan,))
	result = libsakurapy.interpolate_float_yaxis(
		libsakurapy.INTERPOLATION_METHOD_LINEAR, order, nbase, yindata, nchan, zindata, 
		inmask, npos, youtdata, zoutdata, outmask)
	# the result should be [5.25]*nchan
	print 'interpolate_float_yaxis: output {0}'.format(zoutdata)
	del yin, yout, zin, zindata, yindata, youtdata, zoutdata
	# interpolate in X-axis
	xin = [0.,1.]
	xout = [0.25]
	nbase = len(xin)
	npos = len(xout)
	nrow = 3
	zin = [float(6.0), float(5.0)]*nrow
	order = 1
	zindata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (nrow, nbase,))
	zindata[:, 0] = float(6.0)
	zindata[:, 1] = float(5.0)
	xindata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_DOUBLE, (nbase,))
	xindata[:] = xin
	xoutdata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_DOUBLE, (npos,))
	xoutdata[:] = xout
	zoutdata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (nrow, npos,))
	inmask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (nrow, nbase,))
	inmask[:] = True
	outmask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (nrow, npos,))
	result = libsakurapy.interpolate_float_xaxis(libsakurapy.INTERPOLATION_METHOD_LINEAR, order, nbase, xindata, nrow, zindata, inmask, npos, xoutdata, zoutdata, outmask)
	# the result should be [5.75]*nrow
	print 'interpolate_float_xaxis: output {0}'.format(zoutdata)
	
def test_calibration():
	ndata = 7
	yon = [5.0 for i in range(ndata)]
	yoff = [(on - 1.0) for on in yon]
	factor = [float(i) for i in range(ndata)]
	ondata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	ondata[:] = yon
	offdata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	offdata[:] = yoff
	facdata = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	facdata[:] = factor
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	out = libsakurapy.apply_position_switch_calibration(ndata, facdata, ondata, offdata, result)
	# the result should be [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
	print 'calibration: result {0}'.format(out)
	del yon, yoff, factor, ondata, offdata, facdata, result, out
	

def test_convolve1D():
	# Test convolve1d (FFT)
	ndata = 10
	width = 3
	y = [0.]*ndata
	y[5] = 1.0
	peak = ndata / 2
	kernel = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	out = libsakurapy.create_gaussian_kernel_float(peak, width, ndata, kernel)
	print 'convolve1d_fft: kernel {0} (sum {1})'.format(out, out.sum())
	ctx1D = libsakurapy.create_convolve1d_fft_context(ndata, kernel)
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	data[:] = y
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	out = libsakurapy.convolve1d_fft(ctx1D, ndata, data, result)
	print 'convolve1d_fft: input {0} output {1}'.format(data, result)
	del ctx1D, data, result, out
	
	# Test convolve1d (direct)
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	weight = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	mask[:] = True
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	data[:] = y
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	out = libsakurapy.convolve1d(ndata, kernel, ndata, data, mask, result, weight)
	print 'convolve1d: input {0} output {1} weight {2}'.format(data, out, weight)

def test_baseline():
	ndata = 8
	order = 1
	ctxbl = libsakurapy.create_baseline_context(libsakurapy.BASELINE_TYPE_POLYNOMIAL, order, ndata)
	y = [0.5*i for i in range(ndata)]
	m = [True]*ndata
	y[4] += 3.
	m[4] = False
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	data[:] = y
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	mask[:] = m
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))
	final_mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, (ndata,))
	coeff = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_DOUBLE, (order+1,))
	bestfit = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, (ndata,))	
	out = libsakurapy.lsqfit_polynomial(ctxbl, order, ndata, data, mask, 5., 1, order+1, coeff, bestfit, result, final_mask)
	# The result should be [0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0]
	print 'lsqfit_polynomial: input {0}'.format(data)
	print 'lsqfit_polynomial: result {0}'.format(result)
	print 'lsqfit_polynomial: bestfit {0}'.format(bestfit)
	
	# test for bestfit and residual being None
	out2 = libsakurapy.lsqfit_polynomial(ctxbl, order, ndata, data, mask, 5., 1, order+1, coeff, None, None, final_mask)
	print 'lsqfit_polynomial (None version): result {0}'.format(result)
	del ctxbl

def test_complement():
	n = 1024
	dim = (n,)
	data = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, dim)
	mask = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_BOOL, dim)
	data[:] = 10.
	data[5] = numpy.nan
	mask[:] = True
	mask[5] = False
	mask[9] = False
	result = libsakurapy.new_uninitialized_aligned_ndarray(libsakurapy.TYPE_FLOAT, dim)
	libsakurapy.complement_masked_value_float(n, data, mask, result)
	print 'complement_masked_value_float: data {0} mask {1}'.format(data[:10], mask[:10])
	print 'complement_masked_value_float: result {0}'.format(result[:10])
	del data, mask, result

def testAll():
	#test_AB()
	test_stats()
	test_mad()
	test_grid()
	test_logical()
	test_range()
	test_bit()
	test_interpolate()
	test_calibration()
	test_convolve1D()
	test_baseline()
	test_complement()

testAll()
gc.collect(2)

libsakurapy.clean_up()
