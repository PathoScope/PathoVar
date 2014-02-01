import json
import os
import glob

import pathovar
from pathovar.snp_annotation import database_annotator_base, blast_driver
from pathovar.utils.fasta_utils import FastaParser

# Get configuration
database_data = pathovar.get_external_databases_config()['comprehensive_antibiotic_resistance_database']

# Get location of data files
sequence_db_paths = glob.glob(os.path.join(pathovar.INSTALL_DIR, database_data['storage_path'] + '*.fa'))


class CARDNucleotideBlastAnnotator(object):
	def __init__(self, **opts):
		self.opts = opts
		self.verbose = opts.get('verbose', False)
		self.blast_drivers = map(lambda dbf: blast_driver.BlastAnnotationDriver(dbf, "nucleotide"), sequence_db_paths)
		

	@property
	def processes(self):
		return [b.process for b in self.blast_drivers]

	def wait_for_results(self):
		ret_codes =  [b.process.wait() for b in self.blast_drivers]
		if sum(ret_codes) != 0:
			raise blast_driver.BlastDriverException("A CARDNucleotideBlastAnnotator Blast job failed %r" % ret_codes)
		result_files_dict = {self.blast_drivers[ind].db_name:blast_driver.BlastResultsXMLParser(outfile) for ind, outfile in enumerate(self.outfiles)}
		return result_files_dict

	@property
	def outfiles(self):
		return [b.outfile for b in self.blast_drivers]

	def query_with_nucleotides(self, query, **opts):
		for blaster in self.blast_drivers:
			if self.verbose: print("Blasting against %s" % blaster.db_name)
			proc = blaster.blast_with_nucleotides(query, **opts)

def defline_parser(defline):
	pass
