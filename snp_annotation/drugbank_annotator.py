import os
import glob

import pathovar
from pathovar.snp_annotation import blast_driver

class DrugBankNucleotideBlastAnnotator(blast_driver.NucleotideDatabaseBlastAnnotatorBase):
    def __init__(self, storage_path, **opts):
        nucleotide_sequence_db_paths = glob.glob(os.path.join(storage_path, 'nucleotide/*.fasta'))
        blast_driver.NucleotideDatabaseBlastAnnotatorBase.__init__(self, nucleotide_sequence_db_paths, "drugbank", **opts)

class DrugBankProteinBlastAnnotator(blast_driver.ProteinDatabaseBlastAnnotatorBase):
    def __init__(self, storage_path, **opts):
        protein_sequence_db_paths = glob.glob(os.path.join(storage_path, 'protein/*.fasta'))
        blast_driver.ProteinDatabaseBlastAnnotatorBase.__init__(self, protein_sequence_db_paths, "drugbank", **opts)
