import sys
import os
from collections import defaultdict
import argparse

import snp_caller
from snp_caller import snp_utils

argparser = argparse.ArgumentParser(prog="pathovar")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("sam_file", help = "The alignment file to call variants from")
argparser.add_argument('--clean', action = "store_true")
argparser.add_argument("--test", action = "store_true", required = False, help="Enter IPython Interactive Session after execution completes [Development Only]")

target_args = argparser.add_argument_group("Target Organism Selection")
target_args.add_argument("--organism-name", metavar = "ORG-REGEX", action="store", help = "A valid regular expression that matches organism names associated with reference genomes.")
target_args.add_argument("--tax-ids", metavar = "TIs,", default='', action = "store", help = "A comma separated list of NCBI Taxonomy ID numbers for each genome to call against")
target_args.add_argument("--gene-ids", metavar = "GIs,", default='', action = "store", help = "A comma separated list of NCBI Gene ID numbers for each genome to call against -- Experimental")
target_args.add_argument('-r',"--reference-genomes", metavar = "REF", action = "store", required=True, help = "path to a fasta file containing all reference genomes to call against")

snp_caller_args = argparser.add_argument_group("SNP Caller Options")
snp_caller_args.add_argument('-s','--snp-caller', action="store", default = "samtools", choices = ["samtools"], help="Select the SNP Calling Program")
snp_caller_args.add_argument('-b','--snp-caller-binary-location', action="store", default = "", help = "Location of SNP Caller program binaries. Default will search for them on the system path")

snp_anno_args = argparser.add_argument_group("Variant Annotation Options")
snp_anno_args.add_argument('-a', '--annotation-engine', action='store', default = '', choices = ['entrez', 'snpeff', ''], help = "Select the program to annotate variants with")

def main(args):
	snp_caller_driver = None
	opts = {}
	opts['verbose'] = args.verbose
	opts['clean'] = args.clean

	args.tax_ids = args.tax_ids.split(',')
	args.gene_ids = args.gene_ids.split(',')

	if not os.path.exists(args.sam_file): raise IOError("Input .sam File Not Found")
	if args.snp_caller == "samtools":
		from snp_caller import samtools_snp_caller
		snp_caller_driver = samtools_snp_caller.SamtoolsSNPCaller(bin_dir = args.snp_caller_binary_location, opts = opts)

	variant_file = snp_caller_driver.call_snps(args.sam_file, source = args.reference_genomes, org_name = args.organism_name, tax_ids = args.tax_ids, gene_ids = args.gene_ids)

	snp_annotation_driver = None
	if args.annotation_engine == "entrez":	
		from snp_annotation import locate_variant
		snp_annotation_driver = locate_variant.EntrezAnnotationMapper(variant_file, opts)
		snp_annotation_driver.annotate_all_snps()

	if(args.test):
		import IPython
		IPython.embed()

if __name__ == '__main__':
	args = argparser.parse_args()
	if args.verbose: print(args)
	main(args)

