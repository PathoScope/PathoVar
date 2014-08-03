import os
import glob

from pathovar.snp_annotation import blast_driver

class CARDNucleotideBlastAnnotator(blast_driver.NucleotideDatabaseBlastAnnotatorBase):
	def __init__(self, storage_path, bin_dir = '',**opts):
		nucleotide_sequence_db_paths = glob.glob(os.path.join(storage_path, "nucleotide", "*.fa"))
		blast_driver.NucleotideDatabaseBlastAnnotatorBase.__init__(self, nucleotide_sequence_db_paths, "comprehensive_antibiotic_resistance_database", bin_dir,  **opts)

class CARDProteinBlastAnnotator(blast_driver.ProteinDatabaseBlastAnnotatorBase):
	def __init__(self, storage_path, bin_dir = '', **opts):
		protein_sequence_db_paths =  glob.glob(os.path.join(storage_path, "protein", "*.fa"))
		blast_driver.ProteinDatabaseBlastAnnotatorBase.__init__(self, protein_sequence_db_paths, "comprehensive_antibiotic_resistance_database", bin_dir,**opts)



