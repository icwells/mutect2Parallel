'''This script contains functions for ploting values from platypus and mutect2'''

import matplotlib.pyplot as plt

class Columns():

	def __init__(self, header):
		# Defines class for identifying column numbers
		self.ID = None
		self.A = None
		self.B = None
		self.Common = None
		self.Similarity = None
		self.Max = 0
		self.__setColumns__(header)

	def __setColumns__(self, header):
		# Assigns column numbers
		indeces = []
		for idx, i in enumerate(header):
			i = i.strip()
			if i == "ID" or i == "Sample":
				self.ID = idx
				indeces.append(idx)
			elif i == "#PrivateA": 
				self.A = idx
				indeces.append(idx)
			elif i == "#PrivateB":
				self.B = idx
				indeces.append(idx)
			elif i == "#Common":
				self.Common = idx
				indeces.append(idx)
			elif i == "%Similarity" or i == "filtNABcovB_propU":
				self.Similarity = idx
				indeces.append(idx)
			self.Max = max(indeces)

#-----------------------------------------------------------------------------

def getLabels(val):
	# Returns plot lables
	main = "Mutect2/Platypus "
	xlab = "Mutect2: Number of Varaints"
	ylab = "Platypus: Number of Varaints"
	if val == "%":
		main += "Similarities"
		xlab = "Mutect2: Similarity (%)"
		ylab = "Platypus: Similarity (%)"
	elif val == "a":
		main += "Private Variants - A"
	elif val == "b":
		main += "Private Variants - B"
	elif val == "c":
		main += "Variants Common to A and B"
	return main, xlab, ylab

def plotSimilarity(outfile, val, points)
	# Creates scatter plot and calculates regression and correlation
	main, xlab, ylab = getLabels(val)
	fig, ax = plt.subplots()
	ax.scatter(points[0], points[1])
	ax.set_title(main)
	ax.set_xlabel(xlab)
	ax.set_ylabel(ylab)
	fig.savefig(outfile)

