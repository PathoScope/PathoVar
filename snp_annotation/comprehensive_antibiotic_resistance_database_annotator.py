import json
import os
import glob

import pathovar
from pathovar.snp_annotation import database_annotator_base, blast_driver
from pathovar.utils.fasta_utils import FastaParser


database_data = pathovar.get_external_databases_config()['comprehensive_antibiotic_resistance_database']
sequence_db_paths = glob.glob(database_data['storage_path'] + os.sep + '*.fa')




def defline_parser(defline):
    pass
