import os, sys
import vcf

class SNPEffDriver(object):

	def __init__(self, bin_dir, **opts):
		self.bin_dir = bin_dir
		self.opts = opts
		

