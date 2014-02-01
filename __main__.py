

import sys
import os
from time import time
from collections import defaultdict
import argparse

import pathovar
from pathovar import snp_caller
from pathovar import utils

external_database_conf = pathovar.get_external_databases_config()

argparser = argparse.ArgumentParser(prog="pathovar", formatter_class= lambda prog: argparse.HelpFormatter(prog, width=100))
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("sam_file", help = "The alignment file to call variants from [Required]")
argparser.add_argument('--clean', action = "store_true", help="Clean up intermediary files after they're finished being used [default:False]")
argparser.add_argument("--test", action = "store_true", required = False, help="Enter IPython Interactive Session after execution completes [Development Only]")

target_args = argparser.add_argument_group("Target Organism Selection")
target_args.add_argument('-r',"--reference-genomes", metavar = "REF", action = "store", required=True, help = "path to a fasta file containing all reference genomes to call against. [Required]")
target_args.add_argument("--org-names", metavar = "ORG-REGEX", action="store", help = "A valid regular expression that matches organism names associated with reference genomes.[optional]")
target_args.add_argument("--tax-ids", metavar = "TI-REGEX,", action = "store", help = "A valid regular expression that matches NCBI Taxonomy ID numbers for each genome to call against.[optional]")
target_args.add_argument("--gene-ids", metavar = "GI-REGEX,", action = "store", help = "A valid regular expression that matches NCBI Gene ID numbers for each genome to call against.[optional]")


snp_caller_args = argparser.add_argument_group("SNP Caller Options")
snp_caller_args.add_argument('-s','--snp-caller', action="store", default = "samtools", choices = ["samtools"], help="Select the SNP Calling Program.[default:samtools]")
snp_caller_args.add_argument('-b','--snp-caller-binary-location', action="store", default = "", help = "Location of SNP Caller program binaries. Default will search for them on the system path")

snp_anno_args = argparser.add_argument_group("Variant Annotation Options")
snp_anno_args.add_argument('-a', '--annotation-engine', action='store', default = 'entrez', choices = ['entrez', 'snpeff', ''], help = "Select the program to annotate variants with [default:entrez]")
snp_anno_args.add_argument("--no-cache", action= "store_true", default = False, help="Do not use the annotation cache to re-load annotations if they exist to speed up the annotation process.")
snp_anno_args.add_argument("--cache-dir", action = "store", type=str, default='.anno_cache', help="The location to store raw and processed annotation source data. [default='.anno_cache/']")


snp_filt_args = argparser.add_argument_group("Variant Filtering Options")
snp_filt_args.add_argument('--min-depth', type=int, default=5, help="The minimum number of reads that must map to a location to trust a given variant call [default:5]")
snp_filt_args.add_argument('--alt-depth', type=float, default=0.4, help="The minimum ratio of all calls for a locus that must be an alternative allele [default:0.4]")

def main(args):
	opts = {}
	opts['verbose'] = args.verbose
	opts['clean'] = args.clean
	opts['cache_dir'] = args.cache_dir
	opts['no_cache'] = args.no_cache

	if not os.path.exists(args.sam_file): raise IOError("Input .sam File Not Found")
	
	snp_caller_driver = None
	start_clock = time()
	if args.snp_caller == "samtools":
		from snp_caller import samtools_snp_caller
		snp_caller_driver = samtools_snp_caller.SamtoolsSNPCaller(bin_dir = args.snp_caller_binary_location, opts = opts)

	variant_file = snp_caller_driver.call_snps(args.sam_file, source = args.reference_genomes, org_names_reg = args.org_names, tax_ids_reg = args.tax_ids, gene_ids_reg = args.gene_ids)
	consensus_sequences = variant_file + ".cns.fq"

	snp_called_time = time()
	if args.verbose: print('SNP Calling Done (%s sec)' % str(snp_called_time - start_clock))

	snp_annotation_driver = None
	anno_vcf = None

	filter_args = utils.Namespace()
	filter_args.alt_depth = args.alt_depth
	filter_args.min_depth = args.min_depth

	if args.annotation_engine == "entrez":	
		from snp_annotation import locate_variant
		snp_annotation_driver = locate_variant.EntrezAnnotationMapper(variant_file, **dict(filter_args = filter_args, **opts))
		snp_annotation_driver.annotate_all_snps()
		anno_vcf = snp_annotation_driver.write_annotated_vcf()

	snp_annotated_time = time()
	if args.verbose: print('SNP Annotation Done (%s sec)' % str(snp_annotated_time - snp_called_time))

	from pathovar.utils import vcf_utils
	if args.verbose: print("Generating Gene Report.")
	vcf_utils.vcf_to_gene_report(anno_vcf)
	ref_fa = vcf_utils.generate_reference_fasta_for_variants(anno_vcf, snp_annotation_driver.annotation_cache)

	from pathovar.snp_annotation.comprehensive_antibiotic_resistance_database_annotator import CARDNucleotideBlastAnnotator
	card_blast = CARDNucleotideBlastAnnotator()
	card_blast.query_with_nucleotides(ref_fa)


	## Block while annotations run
	

	if(args.test):
		import IPython
		IPython.embed()
	return anno_vcf

if __name__ == '__main__':
	args = argparser.parse_args()
	main(args)






#END