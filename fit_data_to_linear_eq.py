#!/Users/martinmccullagh/anaconda/bin/python

# import libraries
import scipy
from scipy.linalg import *
import sys
import os
import numpy
import math

temperature = 300 # K
kT = numpy.float64(1.9872041E-3*temperature) # kcal/mol
bin_size = numpy.float64(0.01)
thresh = 1e-4

inp_dat_file = sys.argv[1]

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
data = numpy.matrix(data)

rows = data.shape[0]
columns = data.shape[1]

for i in range(columns):
	max_val = numpy.amax(data[:,i])
	min_val = numpy.amin(data[:,i])
	num_bins = int( (max_val-min_val)/bin_size)+1
	prob = numpy.zeros((num_bins,1),dtype=numpy.dtype('f8'))
	x2 = numpy.empty((num_bins,1),dtype=numpy.dtype('f8'))

	avg_x = 0
	for j in range(rows):
		avg_x += data[j,i]
		current_bin = int( (data[j,i] - min_val)/bin_size)
		prob[current_bin,0] += numpy.float64(1.0)

	# finish probability density
	prob /= bin_size*numpy.float64(rows)
	avg_x /= float(rows)

	# compute x2 and take log of probability
	num_nonzero=0
	for j in range(num_bins):
		temp = j* bin_size + min_val - avg_x
		x2[j,0] = numpy.float64(temp*temp)
		if prob[j,0] > thresh :
			# count
			num_nonzero += 1

	# declare arrays that will not contain zeros
	free_energy = numpy.empty((num_nonzero,1),dtype=numpy.dtype('f8'))
	new_prob = numpy.empty((num_nonzero,1),dtype=numpy.dtype('f8'))
	new_x2 = numpy.empty((num_nonzero,1),dtype=numpy.dtype('f8'))
	temp_array = numpy.empty(num_nonzero,dtype=numpy.dtype('f8'))

	# populate nonzero arrays
	num_nonzero = 0
	for j in range(num_bins):
		temp = j* bin_size + min_val - avg_x
		if prob[j,0] > thresh :
			temp_array[num_nonzero] = temp
			free_energy[num_nonzero,0] = numpy.float64(-kT*math.log(prob[j,0]))
			new_prob[num_nonzero,0] = prob[j,0]
			new_x2[num_nonzero,0] = x2[j,0]
			num_nonzero += 1
	
	free_energy -= numpy.amin(free_energy)
	# print some checkpoint stuff
	output_file = "energy_column"+str(i+1)+".dat"
	out = open(output_file,"w")
	for j in range(num_nonzero):
		out.write("%10.5f %10.5f %10.5f %10.5f\n" %(temp_array[j],new_x2[j,0],free_energy[j,0],new_prob[j,0]))
	out.close
	# fit
	k = numpy.linalg.lstsq(new_x2,free_energy)[0]
#	k = scipy.linalg.lstsq(new_x2,free_energy)[0]
	print i+1, k[0][0]
	del prob
	del x2
	del free_energy
	del new_x2
