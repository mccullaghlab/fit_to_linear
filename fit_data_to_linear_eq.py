#!/Users/martinmccullagh/anaconda/bin/python
#
# PROGRAM: fit_data_to_linear_eq.py is a program to take a time series of MD data, compute the probability density, convert that probability density to 
# free energy and then fit that free energy to a harmonic potential.  The output is the force constant (C) for the quadratic equation FE(x) = Cx^2.

# USAGE: fit_data_to_linear_eq.py [config file]

# CONFIG FILE FORMAT:
#	datFile = rmsd.dat
#	temperature = 300
#	thresh = 1e-1
#	binSize = 0.01
#	outFile = rmsd_k_values.dat

# DEPENDENCIES: 
#	scipy
#	numpy


# import libraries
import scipy
from scipy.linalg import *
import sys
import os
import numpy
import math

# Subroutines

# read the configuration file and populate the global variables
def parseConfigFile(cfg_file):
	global inp_dat_file, temperature, bin_size, thresh, out_file
	inp_dat_file = None
	out_file = None
	temperature = 300
	bin_size = 0.01
	thresh = 1E-4
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			# check value
			if option.lower()=='datfile':
				inp_dat_file = value
			elif option.lower()=='outfile':
				out_file = value
			elif option.lower()=='temperature':
				temperature = float(value)
			elif option.lower()=='thresh':
				thresh = float(value)
			elif option.lower()=='binsize':
				binSize = float(value)
			else :
				print "Option:", option, " is not recognized"
	f.close

	if inp_dat_file == None:
		print "No data file is defined.  Please add the following line to your config file:\n"
		print "datFile = [data file name]\n"
		sys.exit("bombing out")
	elif out_file == None:
		print "No output file is defined.  Please add the following line to your config file:\n"
		print "outFile = [output file name]\n"
		sys.exit("bombing out")
	else:
		print "Reading data file:", inp_dat_file
		print "Writing 1/2k values to file:", out_file
		print "Temperature:", temperature
		print "Threshold:", thresh
		print "Bin size:", bin_size

# Main

cfg_file = sys.argv[1]

parseConfigFile(cfg_file)

kT = numpy.float64(1.9872041E-3*temperature) # kcal/mol

# Read data file and populate data matrix (ignoring first column)
f = open(inp_dat_file)
data = []
line_count = 0
counter = 0
for line in f:
	if line_count > 1000:
		data.append([])
		array = line.split()
		for i in range(1,len(array)):
			data[counter].append(numpy.float64(array[i]))
		counter += 1
	line_count += 1
f.close
# transform python list into numpy matrix
data = numpy.matrix(data)
rows = data.shape[0]
columns = data.shape[1]

print "Number of data columns:", columns
print "Number of rows of data:", rows

out = open(out_file, "w")
out.write("#%9s %19s %19s\n" % ("Column", " 1/2k ", "rss"))
out.write("#%9s %19s %19s\n" % (" ", "(kcal/mol/A^2)", " "))
# Loop through each column and compute free energy and fit to a harmonic potential
for i in range(columns):
	# determine domain of data
	max_val = numpy.amax(data[:,i])
	min_val = numpy.amin(data[:,i])
	num_bins = int( (max_val-min_val)/bin_size)+1
	# allocate probability and x2 arrays
	prob = numpy.zeros((num_bins,1),dtype=numpy.dtype('f8'))
	x2 = numpy.empty((num_bins,1),dtype=numpy.dtype('f8'))

	# create histogram of data
	avg_x = 0
	for j in range(rows):
		avg_x += data[j,i]
		current_bin = int( (data[j,i] - min_val)/bin_size)
		prob[current_bin,0] += numpy.float64(1.0)

	# finish probability density
	prob /= bin_size*numpy.float64(rows)
	avg_x /= float(rows)

	# compute x2 and determine how many nonzero (> thresh) values of the probability we have
	num_nonzero=0
	for j in range(num_bins):
		temp = j* bin_size + min_val - avg_x
		x2[j,0] = numpy.float64(temp*temp)
		if prob[j,0] > thresh :
			num_nonzero += 1

	# declare arrays that will not contain zeros
	free_energy = numpy.empty((num_nonzero,1),dtype=numpy.dtype('f8'))
	new_x2 = numpy.empty((num_nonzero,1),dtype=numpy.dtype('f8'))

	# populate nonzero arrays
	num_nonzero = 0
	for j in range(num_bins):
		if prob[j,0] > thresh :
			free_energy[num_nonzero,0] = numpy.float64(-kT*math.log(prob[j,0]))
			new_x2[num_nonzero,0] = x2[j,0]
			num_nonzero += 1
	# Zero the minimum of the free energy	
	free_energy -= numpy.amin(free_energy)
	# fit
	k, rss, rank, s = numpy.linalg.lstsq(new_x2,free_energy)

	out.write("%10d %19.8f %19.8f\n" %(i+1, k[0][0], rss[0]))
	# deallocate arrays
	del prob
	del x2
	del free_energy
	del new_x2

out.close
