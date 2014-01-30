import unittest
import os
import glob

import pathovar
from pathovar.__main__ import main
from pathovar.web import entrez_eutils
from pathovar.snp_caller import samtools_snp_caller
from pathovar.snp_annotation import locate_variant
from pathovar import utils
from pathovar.utils import vcf_utils

global_args = utils.Namespace()
global_args.verbose = True
global_args.clean = True

global_args.org_names = None
global_args.tax_ids = None
global_args.gene_ids = None


global_args.snp_caller = "samtools"
global_args.snp_caller_binary_location = ""

global_args.annotation_engine = "entrez"

global_args.sam_file = "data/updated_CC_287_cDNA_sanger_clean_target0_2.sam"
global_args.reference_genomes = "data/rsv.fa"

global_args.min_depth = 5
global_args.alt_depth = 0.4


class TestFullSamtoolsCallerEntrezAnnotation(unittest.TestCase):
    
    @classmethod
    def setUpClass(self):
        if os.path.basename(os.getcwd()) != 'tests':
            os.chdir('tests')

    def test_step_0_clean_up_files(self):
        print("\nClearing Files")
        try:
            os.remove(*glob.glob('data/*.vcf*'))
            os.remove(*glob.glob('data/*.bam*'))
        except:
            pass


    def test_step_1_call_snps(self):
        print("\nCalling SNPs")
        snp_caller_driver = samtools_snp_caller.SamtoolsSNPCaller(bin_dir = global_args.snp_caller_binary_location, opts = global_args.__dict__)
        snp_caller_driver.call_snps(global_args.sam_file, source = global_args.reference_genomes, org_names_reg = global_args.org_names, tax_ids_reg = global_args.tax_ids, gene_ids_reg = global_args.gene_ids)
        self.assertTrue(os.path.exists(global_args.sam_file[:-3] + 'vcf'), "Calls to samtools|bcftools did not produce a VCF File")
        
    def test_step_2_annotate_snps(self):
        print("\nAnnotating SNPs")
        filter_args = utils.Namespace()
        filter_args.alt_depth = global_args.alt_depth
        filter_args.min_depth = global_args.min_depth
        snp_annotation_driver = locate_variant.EntrezAnnotationMapper(global_args.sam_file[:-3] + 'vcf', **dict(filter_args = filter_args, **global_args.__dict__))
        snp_annotation_driver.annotate_all_snps()
        anno_vcf = snp_annotation_driver.write_annotated_vcf()

        self.assertTrue(os.path.exists(global_args.sam_file[:-3] + 'anno.vcf'), "Calls locate_variant.py did not produce an annotated VCF File")
        vcf_utils.vcf_to_gene_report(anno_vcf)
        self.assertTrue(os.path.exists(global_args.sam_file[:-3] + 'anno.variant_report.tsv'), "Calls to vcf_utils.py did not produce an annotated variant report")
        ref_fa = vcf_utils.generate_reference_fasta_for_variants(anno_vcf, snp_annotation_driver.annotation_cache)
        self.assertTrue(os.path.exists(ref_fa), "Calls to vcf_utils did not produce a reference fasta for sequences with variants")
    



if __name__ == '__main__':
    unittest.main()