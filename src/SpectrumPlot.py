# Author:      Adam Robinson
# Description: This script is designed to plot one or more vdat files 
# containing Raman or PL spectra. It includes cosmic ray removal, peak
# detection and various other functionality.

from datetime        import datetime
from glob            import glob
from scipy.integrate import trapz
from scipy.signal    import find_peaks

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
			self.fields[k.strip()] = float(v.strip())

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

	files = sorted(files)

	# Load the files into a vdat data structure.
	dataset = []
	vdat    = []
	for file in files:
		v = VDATFile(file)
		vdat.append(v)
		dataset.append(np.array(v.data))

	# Calculate the max value across the entire spectrum so we can scale 
	# label locations appropriately.
	max_values = []
	min_values = []
	for data in dataset:
		max_values.append(data[:, 2].max())
		min_values.append(data[:, 2].min())

	max_value = max(max_values)
	min_value = min(min_values)
	_range    = max_value - min_value

	# Remove cosmics if requested.
	if args.correct_cosmics:
		cosmic_centers = []
		cosmic_lefts   = []
		cosmic_rights  = []
		min_prominence = _range / args.prominence_divisor
		for data in dataset:
			probable_cosmics = find_peaks(
				data[:, 2], 
				width=[1, 6], 
				prominence=[min_prominence, 100000]
			)

			cosmic_centers.append(list(probable_cosmics[0]))
			cosmic_lefts.append(list(probable_cosmics[1]['left_ips']))
			cosmic_rights.append(list(probable_cosmics[1]['right_ips']))

	if args.manual_cosmics != []:
		if not args.correct_cosmics:
			cosmic_centers = []
			cosmic_lefts   = []
			cosmic_rights  = []
			for data in dataset:
				cosmic_centers.append([])
				cosmic_lefts.append([])
				cosmic_rights.append([])

		if len(args.manual_cosmics) % 3 != 0:
			print("Manual cosmics must be specified in pairs of three (idx, min, max)")
			exit()
		
		for base_idx in range(len(args.manual_cosmics) // 3):
			idx  = int(args.manual_cosmics[base_idx])
			low  = args.manual_cosmics[base_idx + 1]
			high = args.manual_cosmics[base_idx + 2]

			def index_from_wavenumber(idx, k):
				closest_index = 0
				closest       = 10000
				for i, k0 in enumerate(dataset[idx][:, 1]):
					dist = np.abs(k - k0)
					if dist < closest:
						closest       = dist
						closest_index = i

				return closest_index

			cosmic_centers[idx].append(index_from_wavenumber(idx, (high - low) / 2))
			cosmic_lefts[idx].append(index_from_wavenumber(idx, low))
			cosmic_rights[idx].append(index_from_wavenumber(idx, high))


	# Normalize the data if requested.
	if args.normalize:
		for idx, data in enumerate(dataset):
			integral = trapz(data[:, 2], data[:, 1])
			dataset[idx][:, 2] = dataset[idx][:, 2] / integral


	fig, ax = plt.subplots(1, 1)


	# Offset the data points based on the specified offset.
	if args.offset != 0:
		cumulative_offset = 0
		for idx, data in enumerate(dataset):
			dataset[idx][:, 2] = dataset[idx][:, 2] + cumulative_offset
			cumulative_offset += args.offset

	# Remove the cosmic rays detected in an earlier step (if requested)
	if args.correct_cosmics:
		for idx in range(len(dataset)):
			ds_centers = cosmic_centers[idx]
			ds_rights  = cosmic_rights[idx]
			ds_lefts   = cosmic_lefts[idx]

			for c, l, r in zip(ds_centers, ds_lefts, ds_rights):
				# Find the datapoint that straddles the peak on the
				# left and on the right. We will use these to produce a linear
				# interpolation over the range of the cosmic.
				left_idx    = int(round(l)) - 1
				right_idx   = int(round(r)) + 1
				left_x      = dataset[idx][left_idx,  1]
				right_x     = dataset[idx][right_idx, 1]
				left_point  = dataset[idx][left_idx,  2]
				right_point = dataset[idx][right_idx, 2]
				slope       = (right_point - left_point) / (right_x - left_x)

				# Now we take all of the indices that are greater than left_idx
				# and less than right_idx and replace the data points with an interpolation.
				correction_idx = left_idx + 1
				while correction_idx < right_idx:
					x = dataset[idx][correction_idx, 1]
					dataset[idx][correction_idx, 2] = left_point + (x - left_x) * slope
					correction_idx += 1

	plots = []
	for data in dataset:
		pl, = ax.plot(data[:, 1], data[:, 2])
		plots.append(pl)

	# Label each series based on the specified field in the file header.
	if args.label_fields != []:
		labels = []
		for v in vdat:
			l = []
			for field in args.label_fields:
				l.append(str(v.fields[field]))
			labels.append(" / ".join(l))

		ax.legend(
			plots,
			labels
		)

	# Draw peak labels as requested.
	if args.peaks != []:
		if args.labels == []:
			labels = [str(p) for p in args.peaks]
		else:
			labels = args.labels

		for idx, peak in enumerate(args.peaks):
			ax.axvline(peak, c='green', linestyle='-.')

			# Find the highest value within plus or minus 5 inverse centimeters
			# of this peak and place the label above it.
			highest = 0
			for data in dataset:
				for k, c in zip(data[:, 1], data[:, 2]):
					if k > peak - 5 and k < peak + 5:
						if c > highest:
							highest = c

			ax.text(peak + 2, highest + (_range / 5), labels[idx])

	



	# DEBUG CODE
	# # Plot locations where there are cosmics identified.
	# if args.correct_cosmics:
	# 	for cosmic in cosmics:
	# 		ax.axvline(dataset[0][cosmic, 1], c='red')

	# 	for l, r in zip(lb, rb):
	# 		ax.axvline(dataset[0][int(round(l)), 1], c='green')
	# 		ax.axvline(dataset[0][int(round(r)), 1], c='green')


	#code.interact(local=locals())


	ax.set_xlabel(r"Relative Wavenumber [$cm^{-1}$]")
	if args.normalize:
		ax.set_ylabel(r"Arbitrary Units")
	else:
		ax.set_ylabel(r"Intensity [counts]")
	ax.set_title(args.title)

	plt.show()

	















