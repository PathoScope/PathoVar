import os
import glob

import pathovar
from pathovar.snp_annotation import blast_driver
from pathovar.utils.fasta_utils import FastaParser

# Get configuration
database_data = pathovar.get_external_databases_config()['drugbank']

# Get location of data files
sequence_db_paths = glob.glob(os.path.join(pathovar.INSTALL_DIR, database_data['storage_path'] + '*.fasta'))

class DrugBankNucleotideBlastAnnotator(blast_driver.NucleotideDatabaseBlastAnnotatorBase):
    def __init__(self, **opts):
        blast_driver.NucleotideDatabaseBlastAnnotatorBase.__init__(self, sequence_db_paths, **opts)

# class DrugBankNucleotideBlastAnnotator(object):
#     def __init__(self, **opts):
#         self.opts = opts
#         self.verbose = opts.get('verbose', False)
#         self.blast_drivers = map(lambda dbf: blast_driver.BlastAnnotationDriver(dbf, blast_driver.NUCLEOTIDE), sequence_db_paths)
        

#     @property
#     def processes(self):
#         return [b.process for b in self.blast_drivers]

#     def wait_for_results(self):
#         ret_codes =  [b.process.wait() for b in self.blast_drivers]
#         if sum(ret_codes) != 0:
#             raise blast_driver.BlastDriverException("A DrugBankNucleotideBlastAnnotator Blast job failed %r" % ret_codes)
#         result_files_dict = {self.blast_drivers[ind].db_name:blast_driver.BlastResultsXMLParser(outfile) for ind, outfile in enumerate(self.outfiles)}
#         return result_files_dict

#     @property
#     def outfiles(self):
#         return [b.outfile for b in self.blast_drivers]

#     def query_with_nucleotides(self, query, **opts):
#         for blaster in self.blast_drivers:
#             if self.verbose: print("Blasting against %s" % blaster.db_name)
#             proc = blaster.blast_with_nucleotides(query, **opts)