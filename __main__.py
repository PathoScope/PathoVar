

import sys
import os
from time import time
from collections import defaultdict
import argparse

import pathovar
from pathovar import utils
from pathovar.utils.vcf_utils import EXPOSED_FILTERS

argparser = argparse.ArgumentParser(prog="pathovar", formatter_class= lambda prog: argparse.HelpFormatter(prog, width=100), 
	conflict_handler='resolve')
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("sam_file", help = "The alignment file to call variants from [Required]")
argparser.add_argument('--clean', action = "store_true", help="Clean up intermediary files after they're finished being used [default:False]")
argparser.add_argument("--test", action = "store_true", required = False, help="Enter IPython Interactive Session after execution completes [Development Only]")

target_args = argparser.add_argument_group("Target Organism Selection")
target_args.add_argument('-r',"--reference-genomes", metavar = "REF", action = "store", required=True, help = "path to a fasta file containing all reference genomes to call against. [Required]")
target_args.add_argument("--org-names", metavar = "ORG-REGEX", action="store", help = "A valid regular expression that matches organism names associated with reference genomes.[optional]")
target_args.add_argument("--tax-ids", metavar = "TI-REGEX,", action = "store", help = "A valid regular expression that matches NCBI Taxonomy ID numbers for each genome to call against.[optional]")
target_args.add_argument("--gene-ids", metavar = "GI-REGEX,", action = "store", help = "A valid regular expression that matches NCBI Gene ID numbers for each genome to call against.[optional]")
target_args.add_argument("--keep-all-sequences", action="store_true", default=False, help = "Do NOT discard any sequence in the database that is NOT a complete genome or complete plasmid sequence [optional]")

snp_caller_args = argparser.add_argument_group("SNP Caller Options")
snp_caller_args.add_argument('-s','--snp-caller', action="store", default = "samtools", choices = ["samtools"], help="Select the SNP Calling Program.[default:samtools]")
snp_caller_args.add_argument('-b','--snp-caller-path', action="store", default = "", help = "Location of SNP Caller program binaries. Default will search for them on the system path")

snp_anno_args = argparser.add_argument_group("Variant Annotation Options")
snp_anno_args.add_argument("--cache-dir", action = "store", type=str, default='.anno_cache', help="The location to store raw and processed annotation source data. [default='.anno_cache/']")
snp_anno_args.add_argument('--snpeff-path', default = '', action = 'store', required = False, help = "Path to the snpEff.jar and .config files [default: search system path]")
snp_anno_args.add_argument('--blast-path', default = '', action = 'store', required = False, help = 'Path to the NCBI BLAST+ executables. [default: search system path]')

snp_filt_args = argparser.add_argument_group("Variant Filtering Options")
for filter_type in EXPOSED_FILTERS:
	filter_type.customize_parser(snp_filt_args)

def call_snps(args, **opts):
	snp_caller_driver = None
	if args.snp_caller == "samtools":
		from pathovar.snp_caller import samtools_snp_caller
		snp_caller_driver = samtools_snp_caller.SamtoolsSNPCaller(bin_dir = args.snp_caller_path, **opts)
	variant_file = snp_caller_driver.call_snps(args.sam_file, source = args.reference_genomes, 
		org_names_reg = args.org_names, tax_ids_reg = args.tax_ids, gene_ids_reg = args.gene_ids, 
		keep_all = args.keep_all_sequences)
	consensus_sequences = variant_file + ".cns.fq"
	return variant_file

def main(args):
	opts = {}
	opts['verbose'] = args.verbose
	opts['clean'] = args.clean
	opts['cache_dir'] = args.cache_dir

	if not os.path.exists(args.sam_file): raise IOError("Input .sam File Not Found")
	start_clock = time()

	variant_file = call_snps(args, **opts)
	snp_called_time = time()
	if args.verbose: print('SNP Calling Done (%s sec)' % str(snp_called_time - start_clock))
	variant_locator_driver = None
	anno_vcf = None

	filter_args = utils.Namespace()
	filter_args.alt_depth = args.alt_depth
	filter_args.min_depth = args.min_depth
	filter_args.min_mq = args.min_mq
	filter_args.min_qual = args.min_qual
	filter_args.ref_vcfs = args.ref_vcfs
	filter_args.intersection = args.intersection

	from pathovar.snp_annotation import locate_variant, annotation_report 
	from pathovar.snp_annotation.__main__ import run_annotation_report, find_variant_locations
	from pathovar.web import annotation_manager

	annotation_manager_driver = annotation_manager.EntrezAnnotationManager(**opts)

	anno_vcf = find_variant_locations(variant_file, annotation_manager_driver = annotation_manager_driver,
		**dict(filter_args = filter_args, **opts))

	annotation_report_driver = run_annotation_report(args, anno_vcf, annotation_manager_driver, **opts)

	anno_json = annotation_report_driver.to_json_file()

	from pathovar.visualize.build_html_report import build_report
	build_report(anno_json)

	snp_annotated_time = time()
	if args.verbose: print('SNP Annotation Done (%s sec)' % str(snp_annotated_time - snp_called_time))

	if(args.test):
		import IPython
		IPython.embed()

if __name__ == '__main__':
	args = argparser.parse_args()
	main(args)
