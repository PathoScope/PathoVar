import os
import sys
import argparse
import subprocess
from time import time
from collections import defaultdict

import vcf

import locate_variant



argparser = argparse.ArgumentParser(prog = "snp_annotation")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("vcf_file", help = "The variant call file to annotate")
argparser.add_argument("--test", action = "store_true", required = False)
argparser.add_argument('-d', '--alt-depth', type=float, default=0.4, help="The minimum ratio of all calls for a locus that must be an alternative allele [default:0.4]")

def main():
	args = argparser.parse_args()
	if args.verbose: print(args)
	opts = defaultdict(bool)
	opts['verbose'] = args.verbose
	opts['filter_args'] = args

	if not os.path.exists(args.vcf_file): raise IOError("Input .vcf File Not Found")
	timer = time()
	annotation_mapper = locate_variant.EntrezAnnotationMapper(args.vcf_file, **opts)
	annotation_mapper.annotate_all_snps()

	anno_vcf = annotation_mapper.write_annotated_vcf()
	
	print("Annotation Complete: %s sec" % str(time() - timer))
	if args.test:
		import IPython
		IPython.embed()

	return anno_vcf

if __name__ == '__main__':
	main()