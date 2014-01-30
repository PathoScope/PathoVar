import os
import sys
import argparse
import subprocess
from time import time
from collections import defaultdict

import vcf

import locate_variant

from pathovar import utils

argparser = argparse.ArgumentParser(prog = "snp_annotation")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("vcf_file", help = "The variant call file to annotate")
argparser.add_argument("--test", action = "store_true", required = False)
argparser.add_argument("--no-cache", action= "store_false", default = True, help="Do not use the annotation cache to re-load annotations if they exist to speed up the annotation process.")
argparser.add_argument("--cache-dir", action = "store", type=str, default='.anno_cache', help="The location to store raw and processed annotation source data. [default='.anno_cache/']")
argparser.add_argument('--min-depth', type=int, default=5, help="The minimum number of reads that must map to a location to trust a given variant call [default:5]")
argparser.add_argument('--alt-depth', type=float, default=0.4, help="The MAF threshold, under which variants are ignored [default:0.4]")

def main():
	args = argparser.parse_args()
	if args.verbose: print(args)
	opts = defaultdict(bool)
	opts['verbose'] = args.verbose
	opts['filter_args'] = args

	filter_args = utils.Namespace()
	filter_args.alt_depth = args.alt_depth
	filter_args.min_depth = args.min_depth


	if not os.path.exists(args.vcf_file): raise IOError("Input .vcf File Not Found")
	timer = time()
	annotation_mapper = locate_variant.EntrezAnnotationMapper(args.vcf_file, **opts)
	annotation_mapper.annotate_all_snps()

	anno_vcf = annotation_mapper.write_annotated_vcf()

	from pathovar.utils import vcf_utils
	if args.verbose: print("Generating Gene Report.")
	vcf_utils.vcf_to_gene_report(anno_vcf)

	
	print("Annotation Complete: %s sec" % str(time() - timer))
	if args.test:
		import IPython
		IPython.embed()

	return anno_vcf

if __name__ == '__main__':
	main()