import unittest
import os
import glob

from time import time

import pathovar
from pathovar.__main__ import call_snps, argparser
from pathovar.snp_annotation.__main__ import find_variant_locations, run_annotation_report

from pathovar.web import annotation_manager

from pathovar import utils
from pathovar.utils import vcf_utils


global_args = argparser.parse_args('-r . .'.split())
global_args.verbose = True
global_args.clean = True

global_args.snp_caller = "samtools"
global_args.snp_caller_binary_location = ""

global_args.sam_file = "data/updated_outalign.sam.filt.sam"

global_args.reference_genomes = "data/klebsiella-pneumoniae_ti.fa"
global_args.org_names = None
global_args.gene_ids = None
global_args.tax_ids = "1125630"
global_args.keep_all_sequences = False

global_args.cache_dir = '.anno_cache/'

global_args.min_depth = 5
global_args.alt_depth = 0.4

opts = {}
opts['verbose'] = global_args.verbose
opts['clean'] = global_args.clean
opts['cache_dir'] = global_args.cache_dir

variant_file = None
anno_vcf = None
annotation_manager_driver = None
annotation_report_driver = None
ref_fa = None

class TestFullSamtoolsCallerEntrezAnnotation(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        print(global_args)
        if os.path.basename(os.getcwd()) != 'tests':
            os.chdir('tests')

    def test_step_0_clean_up_files(self):
        print("\nClearing Files")
        try:
            os.remove(*glob.glob('data/*.vcf*'))
            os.remove(*glob.glob('data/*.bam*'))
            os.remove(*glob.glob('data/*.xml*'))
        except:
            pass


    def test_step_1_call_snps(self):
        print("\nCalling SNPs")
        timer = time()
        global variant_file
        variant_file = call_snps(global_args, **opts)
        self.assertTrue(os.path.exists(variant_file), "Calls to call_snps() did not produce a VCF File")
        elapsed = time() - timer
        print("%s ms elapsed" % str(elapsed))
        
    def test_step_2_annotate_snps(self):
        print("\nAnnotating SNPs")
        timer = time()
        global annotation_manager_driver
        annotation_manager_driver = annotation_manager.EntrezAnnotationManager(**global_args.__dict__)
        global anno_vcf
        anno_vcf = find_variant_locations(variant_file, annotation_manager_driver, **dict(filter_args = global_args, **opts))
        self.assertTrue(os.path.exists(anno_vcf), "Calls find_variant_locations() did not produce an annotated VCF File")
        elapsed = time() - timer
        print("%s ms elapsed" % str(elapsed))
    
    def test_step_3_external_database_search(self):
        timer = time()
        global annotation_report_driver
        annotation_report_driver = run_annotation_report(global_args, anno_vcf, annotation_manager_driver, **opts)
        result_json = annotation_report_driver.to_json_file()
        self.assertTrue(os.path.exists(result_json), "Calls to run_annotation_report() did not produce an a JSON File")
        elapsed = time() - timer
        print("%s ms elapsed" % str(elapsed))

if __name__ == '__main__':
    unittest.main()
