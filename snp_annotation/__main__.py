import os
import sys
import argparse
import subprocess
from time import time

import vcf

import pathovar
from pathovar import utils
from pathovar.snp_annotation import locate_variant

argparser = argparse.ArgumentParser(prog = "snp_annotation")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("vcf_file", help = "The variant call file to annotate")
argparser.add_argument("--test", action = "store_true", required = False)
argparser.add_argument("--no-cache", action= "store_true", default = False, help="Do not use the annotation cache to re-load annotations if they exist to speed up the annotation process.")
argparser.add_argument("--cache-dir", action = "store", type=str, default='.anno_cache', help="The location to store raw and processed annotation source data. [default='.anno_cache/']")
argparser.add_argument('--min-depth', type=int, default=5, help="The minimum number of reads that must map to a location to trust a given variant call [default:5]")
argparser.add_argument('--alt-depth', type=float, default=0.4, help="The MAF threshold, under which variants are ignored [default:0.4]")

def main():
	args = argparser.parse_args()
	if args.verbose: print(args)
	opts = dict()
	opts['verbose'] = args.verbose
	opts['cache_dir'] = args.cache_dir
	opts['no_cache'] = args.no_cache

	filter_args = utils.Namespace()
	filter_args.alt_depth = args.alt_depth
	filter_args.min_depth = args.min_depth

	opts['filter_args'] = filter_args

	if not os.path.exists(args.vcf_file): raise IOError("Input .vcf File Not Found")
	timer = time()
	annotation_mapper = locate_variant.EntrezAnnotationMapper(args.vcf_file, **opts)
	annotation_mapper.annotate_all_snps()

	anno_vcf = annotation_mapper.write_annotated_vcf()

	from pathovar.utils import vcf_utils
	if args.verbose: print("Generating Gene Report.")
	vcf_utils.vcf_to_gene_report(anno_vcf)
	ref_fa = vcf_utils.generate_reference_fasta_for_variants(anno_vcf, annotation_mapper.annotation_cache)

	## Load internal configuration file
	external_database_conf = pathovar.get_external_databases_config()
	enabled_databases = [database_name for database_name, database_conf in external_database_conf.items() if database_conf['enabled']]
	external_database_results = {}
	# OPT-IN DATABASES
	if "comprehensive_antibiotic_resistance_database" in enabled_databases:
		from pathovar.snp_annotation.comprehensive_antibiotic_resistance_database_annotator import CARDNucleotideBlastAnnotator
		card_blast = CARDNucleotideBlastAnnotator()
		card_blast.query_with_nucleotides(ref_fa)

		# Block while annotations run
		external_database_results["comprehensive_antibiotic_resistance_database"] = card_blast.wait_for_results()

	if "drugbank" in enabled_databases:
		from pathovar.snp_annotation.drugbank_annotator import DrugBankNucleotideBlastAnnotator
		drugbank_blast = DrugBankNucleotideBlastAnnotator()
		drugbank_blast.query_with_nucleotides(ref_fa)

		# Block while annotations run
		external_database_results["drugbank"] = drugbank_blast.wait_for_results()

	if(args.test):
		import IPython
		IPython.embed()
if __name__ == '__main__':
	main()
