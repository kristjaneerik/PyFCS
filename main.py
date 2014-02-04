"""
The handling of FCS files is based on the package FlowPy
(http://flowpy.wikidot.com) by Revanth Sai Kumar, Tejas Mehta, Biplab Bose. It's
been modified to use numpy and sees a ~4.5x speed improvement.

The current version is bound to be full of bugs but is nice for some quick-and-dirty
plotting of FCS files.  Hence this code is provided as-is, under the GPL v2 license
(see LICENSE).
"""

verbose = True
# in the filename formats, %s replaces the input fcs filename (without .fcs)
system_parameters_filename_format = "%s-parameters.dat"
raw_output_filename_format = "%s-raw.csv"
gated_output_filename_format = "%s-gated.csv"
logs = [ 1, 2, 3, 4, 5, 6, 7, 8, 9 ]

from os.path import basename, splitext, isfile
from os import listdir
from numpy import reshape, logspace, insert, array, linspace, log10, concatenate, histogram, searchsorted
from struct import unpack
import matplotlib.pyplot as plt
from operator import mul
from matplotlib.path import Path

def opportunistic_number(potential_number):
	"""
	Formats a string into a number (int, float) if it is one. If it isn't just
	returns the string.
	"""
	try:
		return int(potential_number)
	except ValueError: # wasn't an int
		try:
			return float(potential_number)
		except ValueError: # wasn't a float either
			return potential_number

def parse_FCS_file(input_filename, dump_params = False):
	"""
	Reads a .fcs file. Returns a tuple
		(input_name, system_parameters, parameter_names, datamatrix)
	where input_name is a string holding the FCS filename without the .fcs part,
	system_parameters is a dictionary mapping parameter names to values,
	parameter_names is a list of all the parameters found in the .fcs file
	datamatrix is a matrix of all events and the recorded values.
	If dump_params is True, then dump the parameters into a file named by
	system_parameters_filename_format.
	"""
	if verbose:
		print "Analyzing", input_filename
	
	f = open(input_filename, 'rb')
	input_name = splitext(basename(input_filename))[0]
	
	# Read FCS file parameters (where text starts and ends, data starts and ends etc.)
	fileHeaderver = f.read(10) # reads 10 bytes
	fileHeadertextStart = int(f.read(8))
	fileHeadertextEnd = int(f.read(8))
	fileHeaderdataStart = int(f.read(8))
	fileHeaderdataEnd = int(f.read(8))
	fileHeaderanalysisStart =int(f.read(8))
	fileHeaderanalysisEnd = int(f.read(8))

	f.seek(fileHeadertextStart,0) # move pointer back to beginning
	
	# read header, cleans up chr(0)-s
	fileHeadertextData = f.read(fileHeadertextEnd - fileHeadertextStart + 1).replace(chr(0), '')
	# get list of parameter names, values
	textCell = fileHeadertextData.split(fileHeadertextData[0]) # the first character is chr(12), which is essentially a lot of whitespace that separates parameter names and values that come in pairs
	del textCell[0] # first item in the list is a blank, remove it
	
	# save system parameters in a file and a dictionary for later use
	system_parameters = {}
	if dump_params:
		with open(system_parameters_filename_format % input_name, "w") as fparams:
			for i in xrange(0, len(textCell) - 1, 2):
				fparams.write("%s\t%s\n" % (textCell[i], textCell[i + 1]))
	for i in xrange(0, len(textCell) - 1, 2):
		system_parameters[textCell[i]] = opportunistic_number(textCell[i + 1])
	
	num_events = system_parameters['$TOT'] # number of events recorded
	num_params = system_parameters['$PAR'] # number of data parameters (number of datapoints per event, e.g. things like FSC-A, FITC)
	datatype = system_parameters['$DATATYPE']
	byteord = system_parameters['$BYTEORD']
	bits = system_parameters['$P1B']
	
	f.seek(fileHeaderdataStart, 0) # go to the part of the file where the data starts
	num = fileHeaderdataEnd - fileHeaderdataStart + 1
	data = f.read(num) # read the data
	
	datalist = [None] * (num_events * num_params)
	
	j = 0
	
	# TODO: there are other dataypes, byteords that FlowPy supports -- should
	# find some example files to test; until then, this works with files produced
	# by the LSR Fortessa
	if datatype == 'F':
		if byteord == '1,2,3,4':
			for i in range(0, num / 4):
				datalist[i] = unpack('f', data[4 * i : 4 * i + 4])[0]
		if byteord == '4,3,2,1':
			for i in range(0, num / 4):
				datalist[i] = unpack('f', data[4 * i : 4 * i + 4][::-1])[0]
		else:
			print "Unsupported byte order in FCS file", byteord
	else:
		print "Unsupported datatype in FCS file:", datatype
	
	datamatrix = reshape(datalist, (num_events, num_params))
	
	parameter_names = [ system_parameters["$P%dN" % (i + 1) ] for i in xrange(num_params) ] # things like FSC-A, FITC; they are 0 indexed in the FCS file, ($P1N is the name of the first parameter)
	
	kf = 0
	datal = [1] * (num_events * num_params)

	#logarithmic data # TODO: understand; does the LSR produce data like this?
	loge = [0] * (num_params) # bit list to determine whether a parameter is logarithmic
	for i in range(1, num_params + 1):
		decade = int(system_parameters['$P%dE' % i][0])
		if decade != 0:
			print "I don't understand this quite yet!"
			Range = system_parameters['$P%dR' % i]
			i = i -1
			m = 0
			loge[i] = 1
			for j in range(kf * num_events, num_events + kf * num_events):
				datal[j] = float(10**(float((datam2[0][m][i])*float(decade)/float(Range))))
				m = m + 1
			kf = kf + 1
	if kf != 0: # TODO: understand
		datal = reshape(datal, num_events, num_params)
	j = 0
	
	#Parameters which have log values # TODO: understand
	Paramsl = [None] * (kf)
	for i in range(0, num_params):
		if loge[i] == 1: # this parameter is logarithmic
			Paramsl[j] = parameter_names[i] + 'l'
			j = j + 1

	return input_name, system_parameters, parameter_names, datamatrix

def smooth(hist, level):
	smoothed = list(hist[0:level])
	for i in xrange(level, len(hist)):
		if True:# or i < lin_num_bins:
			smoothed.append(sum(hist[i - level : i + level]) / (2 * level + 1.0))
		else:
			smoothed.append(reduce(mul, hist[i - level : i + level], 1) ** (1 / (2 * level + 1.0)))
	return smoothed

class FCS_Data():
	def __init__(self, filename, write_files = False):
		self.input_name, self.system_parameters, self.parameter_names, \
			self.datamatrix = parse_FCS_file(filename, write_files)
		self.filename = filename
		if verbose:
			print "FCS file contains the parameters:", self.parameter_names
		self.param_index = dict([ (self.parameter_names[i], i) for i in xrange(len(self.parameter_names)) ])
		self.is_gated = False
		self.num_params = len(self.parameter_names)
		
		if write_files:
			self.to_csv()
	
	def set_gating(self, gate_file):
		"""
		Retrieves the gating from an XML file (gate_file) that can be obtained
		by copying a gate in FlowJo and pasting it into a text file.
		"""
		with open(gate_file, 'r') as f:
			gate_text = f.read()
			for bad in [ "gating:", "data-type:" ]:
				gate_text = gate_text.replace(bad, "")
		from xml.dom import minidom
		xml_gate = minidom.parseString(gate_text)
		gates = xml_gate.getElementsByTagName('PolygonGate')
		FSCSSC_A = []
		FSC_HW = []
		SSC_HW = []
		# NOTE: assumes that pairs are *-H, *-W!
		for gate in gates:
			axes = []
			for d in gate.getElementsByTagName('dimension'):
				axes.append(d.getElementsByTagName('parameter')[0].attributes['name'].value)
			
			for vertex in gate.getElementsByTagName('vertex'):
				coordinates = []
				for coordinate in vertex.getElementsByTagName('coordinate'):
					coordinates.append(float(coordinate.attributes['value'].value))
				if axes == [ 'FSC-A', 'SSC-A' ]:
					FSCSSC_A.append(coordinates)
				elif axes == [ 'FSC-H', 'FSC-W' ]:
					FSC_HW.append(coordinates)
				elif axes == [ 'SSC-H', 'SSC-W' ]:
					SSC_HW.append(coordinates)
				else:
					print "Unknown axes for setting gate:", axes
		
		FSCSSC_path = Path(FSCSSC_A)
		FSC_path = Path(FSC_HW)
		SSC_path = Path(SSC_HW)
		"""
		# this code plots the datapoints in the FSC-H/W plane as well as the gating
		fig, ax = plt.subplots()
		ax.set_xscale('log')
		plt.plot(list(self.datamatrix[:,self.param_index['FSC-H']]), list(self.datamatrix[:,self.param_index['FSC-W']]), 'bo')
		from matplotlib import patches
		ax.add_patch(patches.PathPatch(SSC_path, facecolor = 'green'))
		"""
		#return
		
		p1 = self.datamatrix[:,[ self.param_index['FSC-A'], self.param_index['SSC-A'] ]]
		p2 = self.datamatrix[:,[ self.param_index['FSC-H'], self.param_index['FSC-W'] ]]
		p3 = self.datamatrix[:,[ self.param_index['SSC-H'], self.param_index['SSC-W'] ]]
		first = FSCSSC_path.contains_points(p1)
		second = FSC_path.contains_points(p2)
		third = SSC_path.contains_points(p3)
		combined = first & second & third
		gated_rows = []
		gated_data = []
		for i in xrange(len(first)):
			if combined[i] == True:
				gated_rows.append(i)
				gated_data.append(self.datamatrix[i:i+1][0])
		self.gated_rows = array(gated_rows)
		self.gated_data = array(gated_data)
		self.is_gated = True
	
	def histogram(self, parameter):
		"""
		Plots a histogram. Uses gated data if a gate has been set, full data
		otherwise.
		"""
		max_linear = 10
		lin_num_bins = 10
		log_num_bins = 500
		if not self.is_gated:
			print "You haven't gated the population!"
			data = self.datamatrix[:,self.param_index[parameter]]
		else:
			data = self.gated_data[:,self.param_index[parameter]]
		
		fig, ax = plt.subplots()
		ax.set_xscale('symlog', linthreshx = max_linear, subsx = logs )
		linear_bins = linspace(0, max_linear, lin_num_bins, endpoint = False)
		log_bins = logspace(log10(max_linear), 6, log_num_bins)
		bins = concatenate((linear_bins, log_bins))
		#n, bins, patches = plt.hist(data, bins = bins, edgecolor = 'none')
		
		hist, edges = histogram(data, bins = bins)
		
		
		plt.plot(edges[:-1], hist, alpha = 0.1)
		#plt.plot(edges[:-1], smooth(hist, 1))
		#plt.plot(edges[:-1], smooth(hist, 2))
		plt.plot(edges[:-1], smooth(hist, 10))
		#plt.plot(edges[:-1], smooth(hist, 20))
		
		ax.set_xlabel(parameter)
		ax.set_ylabel("Count")
		plt.show()
	
	def mean_ratio(self, numerator, denominator):
		"""
		Finds ratio of numerator / denominator where these are the names of
		parameters that were recorded for events (e.g. FITC-A, Pacific Blue-A)
		"""
		if not self.is_gated:
			print "You haven't gated the population!"
			numerators = self.datamatrix[:,self.param_index[numerator]]
			denominators = self.datamatrix[:,self.param_index[denominator]]
		else:
			numerators = self.gated_data[:,self.param_index[numerator]]
			denominators = self.gated_data[:,self.param_index[denominator]]
		# plot all data
		
		#fig, ax = plt.subplots()
		#ax.set_xscale('symlog', linthreshx = 100, subsx = logs )
		#ax.set_yscale('symlog', linthreshy = 100, subsy = logs )
		#plt.plot(list(numerators), list(denominators), '.', markersize = 1)
		#plt.show()
		
		# figure out histogram for the denominator
		max_linear = 10
		lin_num_bins = 10
		log_num_bins = 500
		smoothing_level = 10
		linear_bins = linspace(0, max_linear, lin_num_bins, endpoint = False)
		log_bins = logspace(log10(max_linear), 6, log_num_bins)
		bins = concatenate((linear_bins, log_bins))
		hist, edges = histogram(denominators, bins = bins)
		smoothed_hist = smooth(hist, smoothing_level)
		
		numerator_sums = [ 0.0 ] * len(edges)
		numerator_counts = [ 0.0 ] * len(edges)
		
		for i in xrange(len(numerators)):
			den = denominators[i]
			den_bin = searchsorted(edges, den)
			#print numerators[i], den_bin, den, edges[den_bin - 1], edges[den_bin], edges[den_bin + 1]
			numerator_sums[den_bin] += numerators[i]
			numerator_counts[den_bin] += 1
		
		numerator_averages = array(numerator_sums) / array(numerator_counts)
		
		#fig, ax = plt.subplots()
		#ax.set_xscale('symlog', linthreshx = 100, subsx = logs )
		#ax.set_yscale('symlog', linthreshy = 100, subsy = logs )
		normalized = numerator_averages / edges
		smooth_normalized = smooth(normalized, 15)
		#print smooth_normalized # many nans
		#plt.plot(edges, smooth_normalized, 'r')
		#plt.show()
		return edges, smooth_normalized
	
	def moving_avg_relative(self, numerator, denominator):
		"""
		Finds ratio of numerator / denominator where these are the names of
		parameters that were recorded for events (e.g. FITC-A, Pacific Blue-A).
		Uses a moving average, averaging values 0.85 ... 1.15 times the current
		value.
		"""
		from operator import itemgetter
		
		if not self.is_gated:
			print "You haven't gated the population!"
			numerators = self.datamatrix[:,self.param_index[numerator]]
			denominators = self.datamatrix[:,self.param_index[denominator]]
		else:
			numerators = self.gated_data[:,self.param_index[numerator]]
			denominators = self.gated_data[:,self.param_index[denominator]]
		
		# create pairs of numerator, denominator and sort them by the denominator value
		pairs = sorted(zip(numerators, denominators), key = itemgetter(1))
		
		N = len(pairs)
		avg_nums = []
		avg_denoms = []
		for i in xrange(N):
			# TODO inefficient but meh
			num = pairs[i][0]
			s_num = num
			den = pairs[i][1]
			s_den = den
			j = i + 1
			top = 1.15 * den
			while j < N and pairs[j][1] <= top:
				s_num += pairs[j][0]
				s_den += pairs[j][1]
				j += 1
			k = i - 1
			bottom = 0.85 * den
			while k >= 0 and pairs[k][1] >= bottom:
				s_num += pairs[k][0]
				s_den += pairs[k][1]
				k -= 1
			count = float(j - k + 1)
			avg_nums.append(s_num / count)
			avg_denoms.append(s_den / count)
			if verbose and i % 100 == 0:
				print i, i / float(N)
		
		avg_nums = array(avg_nums)
		avg_denoms = array(avg_denoms)
		
		return avg_nums, avg_denoms
	
	def to_csv(self):
		"""
		Writes the raw data into a csv file. If a gating is set, output only
		the gated population.
		"""
		if not self.is_gated:
			data = self.datamatrix
			filename_format = raw_output_filename_format
		else:
			data = self.gated_data
			filename_format = gated_output_filename_format
		
		FILE = open(filename_format % self.input_name, "w")
		FILE.write('\t')
		FILE.write("\t".join(self.parameter_names)) # write parameter names
		# TODO: so-called log parameters
		#for t in range(0, self.num_params):
		#	if loge[t] == 1: #log parameters
		#		FILE.write("%s (log)\t" % parameter_names[t])
		FILE.write('\n')
		
		row_format = "\t%d" * self.num_params
		rows = []
		num_events = len(data)
		for i in xrange(num_events):
			rows.append(row_format % tuple(data[i]))
			# TODO: log values (things in datal)
		FILE.write("\n".join(rows)) # 0.1s faster than writing every line separately for 16k events
		FILE.close()
		"""
		#Printing Data values in the file # TODO: understand
		for i in range(0, num_events):
			s1 = str('\t',)
			FILE.write(s1)
			for j in range(0, num_params):
				FILE.write("%d\t" % datamatrix[i][j])
				k = 0
				if j == num_params-1:
					for t in range(0, num_params): #log values
						if loge[t] == 1:
							k = k + 1
							s = str(datal[i][k-1],)
							FILE.write(s)
							FILE.write('\t')
			s2 = '\n'
			FILE.write(s2)
		"""

def plot_pairs(directory, gate_file, pairs, name_map, c1, c2):
	for pair in pairs:
		sample1 = FCS_Data(directory + "/" + pair[0])
		sample1.set_gating(directory + "/" + gate_file)
		y1, x1 = sample1.moving_avg(c1, c2)
		
		sample2 = FCS_Data(directory + "/" + pair[1])
		sample2.set_gating(directory + "/" + gate_file)
		y2, x2 = sample2.moving_avg(c1, c2)
		
		fig, ax = plt.subplots()
		plt.plot(x1, y1, '.', label = name_map[pair[0]])
		plt.plot(x2, y2, '.', label = name_map[pair[1]])
		ax.set_ylabel(c1)
		ax.set_xlabel(c2)
		ax.set_xscale('symlog', linthreshx = 100, subsx = logs)
		ax.set_yscale('symlog', linthreshy = 100, subsy = logs)
		
		plt.legend(loc = 'upper left')
		#plt.show()
		plt.savefig("%s/%s %s.png" % (directory, pair[0], pair[1]), dpi = 200)

if __name__ == "__main__":
	from sys import argv, exit
	
	"""
	# example 1: load FCS file, apply a gate and save to CSV file
	Sample6 = FCS_Data("Samples_006.fcs")
	Sample6.set_gating("gate.xml")
	Sample6.to_csv()
	"""
	
	"""
	# example 2: convert all FCS files in a directory to CSV files; uses the gate.xml file in the folder to apply a gating
	folder = "." # at least on Linux and Mac, . means the current folder
	# folder = "./data" # another example of setting the folder
	all_files = [ f for f in listdir(folder) if isfile(f) ]
	fcs_files = [ f for f in all_files if splitext(basename(f))[1] in [ ".fcs", ".FCS" ] ]
	for f in fcs_files:
		D = FCS_Data(f)
		D.set_gating(folder + "/gate.xml")
		D.to_csv()
	"""
	
	# example 3: TODO
	# TODO
