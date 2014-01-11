import os
import sys
import argparse
import subprocess
from collections import defaultdict
import locate_variant

argparser = argparse.ArgumentParser(prog = "snp_annotation")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("vcf_file", help = "The variant call file to annotate")
argparser.add_argument("--test", action = "store_true", required = False)

def main():
	args = argparser.parse_args()
	if args.verbose: print(args)
	opts = defaultdict(bool)
	opts['verbose'] = args.verbose

	if not os.path.exists(args.vcf_file): raise IOError("Input .vcf File Not Found")
	
	annotation_mapper = locate_variant.EntrezAnnotationMapper(args.vcf_file, opts)
	annotation_mapper.annotate_all_snps()

	annotation_mapper.write_annotated_vcf()

	if args.test:
		import IPython
		IPython.embed()

if __name__ == '__main__':
	main()