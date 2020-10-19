# Author:      Adam Robinson
# Description: This script is designed to plot one or more vdat files 
# containing Raman or PL spectra. It includes cosmic ray removal, peak
# detection and various other functionality.

from datetime import datetime
from glob     import glob

import numpy             as np
import matplotlib.pyplot as plt

import argparse
import code
import json

# Parses a vdat file and extracts the timestamp, header values and
# data. The data is stored in VDATfile.data, where the columns are
#
# 1) Wavelength (nm)
# 2) Wavenumber (cm^-1)
# 3) Intensity (counts)
class VDATFile:
	def __init__(self, path):
		with open(path, 'r') as file:
			self.raw_text = file.read()

		self.raw_text.replace('\r', '')
		self.lines = self.raw_text.split("\n")
		self.lines = [l for l in self.lines if l.strip() != '']

		self.parseHeader()
		self.parseData()

	def parseHeader(self):
		# First, read the date and time.
		self.datetime = datetime.strptime(self.lines[0], '%m/%d/%Y %I:%M %p')

		# The second line should be the fields.
		raw_fields  = self.lines[1].split("\t")
		self.fields = {} 

		for f in raw_fields:
			k, v = f.split(':')
			self.fields[k.strip().lower()] = float(v.strip())

		self.comments = self.lines[2].split(':')[1].strip()

	def parseData(self):
		self.data = [list(map(float, line.split('\t'))) for line in self.lines[5:]]


# Process the command line arguments supplied to the program. These will be 
# in the json file named "SpectrumPlot.py"
def preprocess(args_specification):
	parser = argparse.ArgumentParser(description=args_specification['description'])

	types = {'str': str, 'int': int, 'float': float}

	for argument in args_specification['arguments']:
		spec = argument['spec']
		if 'type' in spec:
			spec['type'] = types[spec['type']]
		parser.add_argument(
			*argument['names'], 
			**spec
		)

	args = parser.parse_args()

	return args

if __name__ == '__main__':
	with open("SpectrumPlot.json", 'r') as file:
		args_specification = json.loads(file.read())

	args = preprocess(args_specification)

	# Load all of the files specified by the user.
	files = []
	for path in args.inputs:
		files.extend(glob(path))

	fig, ax = plt.subplots(1, 1)

	for file in files:
		data = np.array(VDATFile(file).data)
		ax.plot(data[:, 1], data[:, 2])

	ax.set_xlabel(r"Relative Wavenumber [$cm^{-1}$]")
	ax.set_ylabel(r"Intensity [counts]")
	ax.set_title("Spectrum")

	plt.show()

	















