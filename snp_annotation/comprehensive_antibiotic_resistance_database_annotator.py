import json
import os
import glob

import pathovar
from pathovar.snp_annotation import database_annotator_base, blast_driver
from pathovar.utils.fasta_utils import FastaParser

# Get configuration
database_data = pathovar.get_external_databases_config()['comprehensive_antibiotic_resistance_database']

# Get location of data files
nucleotide_sequence_db_paths = glob.glob(os.path.join(pathovar.INSTALL_DIR, database_data['storage_path'] + 'nucleotide/*.fa'))
protein_sequence_db_paths = glob.glob(os.path.join(pathovar.INSTALL_DIR, database_data['storage_path'] + 'protein/*.fa'))



class CARDNucleotideBlastAnnotator(blast_driver.NucleotideDatabaseBlastAnnotatorBase):
	def __init__(self, storage_path, **opts):
		nucleotide_sequence_db_paths = glob.glob(os.path.join(storage_path, "nucleotide", "*.fa"))
		blast_driver.NucleotideDatabaseBlastAnnotatorBase.__init__(self, nucleotide_sequence_db_paths, "comprehensive_antibiotic_resistance_database", **opts)

class CARDProteinBlastAnnotator(blast_driver.ProteinDatabaseBlastAnnotatorBase):
	def __init__(self, storage_path, **opts):
		protein_sequence_db_paths =  glob.glob(os.path.join(storage_path, "protein", "*.fa"))
		blast_driver.ProteinDatabaseBlastAnnotatorBase.__init__(self, protein_sequence_db_paths, "comprehensive_antibiotic_resistance_database", **opts)

def defline_parser(defline):
	pass

if __name__ == '__main__':
	CARDNucleotideBlastAnnotator()
	CARDProteinBlastAnnotator()
