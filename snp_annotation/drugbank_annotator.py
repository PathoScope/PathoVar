import os
import glob

import pathovar
from pathovar.snp_annotation import blast_driver
from pathovar.utils.fasta_utils import FastaParser

# Get configuration
database_data = pathovar.get_external_databases_config()['drugbank']

# Get location of data files
nucleotide_sequence_db_paths = glob.glob(os.path.join(pathovar.INSTALL_DIR, database_data['storage_path'] + 'nucleotide/*.fasta'))
protein_sequence_db_paths = glob.glob(os.path.join(pathovar.INSTALL_DIR, database_data['storage_path'] + 'protein/*.fasta'))

class DrugBankNucleotideBlastAnnotator(blast_driver.NucleotideDatabaseBlastAnnotatorBase):
    def __init__(self, **opts):
        blast_driver.NucleotideDatabaseBlastAnnotatorBase.__init__(self, nucleotide_sequence_db_paths, "drugbank", **opts)

class DrugBankProteinBlastAnnotator(blast_driver.ProteinDatabaseBlastAnnotatorBase):
    def __init__(self, **opts):
        blast_driver.ProteinDatabaseBlastAnnotatorBase.__init__(self, protein_sequence_db_paths, "drugbank", **opts)
