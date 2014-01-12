import unittest
import os
import glob

import pathovar
from pathovar.__main__ import main
from pathovar.web import entrez_eutils
from pathovar.snp_caller import samtools_snp_caller
from pathovar.snp_caller import snp_utils

global_args = snp_utils.Namespace()
global_args.verbose = True
global_args.clean = True
#global_args.test = True

global_args.org_names = None
global_args.tax_ids = None
global_args.gene_ids = None


global_args.snp_caller = "samtools"
global_args.snp_caller_binary_location = ""

global_args.annotation_engine = "entrez"

global_args.sam_file = "data/updated_CC_287_cDNA_sanger_clean_target0_2.sam"
global_args.reference_genomes = "data/rsv.fa"

class TestFullSamtoolsCallerEntrezAnnotation(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        try:
            os.remove(*glob.glob('data/*.vcf*'))
        except:
            pass

    def test_step_1_call_snps(self):
        print("\nCalling SNPs")
        snp_caller_driver = samtools_snp_caller.SamtoolsSNPCaller(bin_dir = global_args.snp_caller_binary_location, opts = global_args.__dict__)
        snp_caller_driver.call_snps(global_args.sam_file, source = global_args.reference_genomes, org_names_reg = global_args.org_names, tax_ids_reg = global_args.tax_ids, gene_ids_reg = global_args.gene_ids)
        self.assertTrue(os.path.exists(global_args.sam_file[:-3] + 'vcf'), "Calls to samtools|bcftools did not produce a VCF File")
        
    def test_step_2_annotate_snps(self):
        print("\nAnnotating SNPs")

if __name__ == '__main__':
    unittest.main()