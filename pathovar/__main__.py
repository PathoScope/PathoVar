

import sys
import os
from time import time
from collections import defaultdict
import argparse

import pathovar
from pathovar import utils
from pathovar.utils import vcf_utils
from pathovar.utils.config import DEFAULT_CONFIG_PATH, get_config, load_param

from pathovar.utils.vcf_utils import EXPOSED_FILTERS
from pathovar.utils.vcf_utils import FilterByAltCallDepth
from pathovar.utils.vcf_utils import FilterByReadDepth
from pathovar.utils.vcf_utils import FilterByMappingQuality
from pathovar.utils.vcf_utils import FilterByCallQuality
from pathovar.utils.vcf_utils import FilterByComparisonVCF

argparser = argparse.ArgumentParser(prog="pathovar", formatter_class= lambda prog: argparse.HelpFormatter(prog, width=100), 
	conflict_handler='resolve')
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("sam_file", help = "The alignment file to call variants from [Required]")
argparser.add_argument('--clean', action = "store_true", help="Clean up intermediary files after they're finished being used [default:False]")
argparser.add_argument("-c", "--config", action = "store", default = None, help = "Path to the pathovar.conf.json configuration file to use [default:System Wide]")
argparser.add_argument("--test", action = "store_true", required = False, help="Enter IPython Interactive Session after execution completes [Development Only]")

target_args = argparser.add_argument_group("Target Organism Selection")
target_args.add_argument('-r',"--reference-genomes", metavar = "REF", action = "store", required=True, help = "path to a fasta file containing all reference genomes to call against. [Required]")
target_args.add_argument("--org-names", metavar = "ORG-REGEX", action="store", help = "A valid regular expression that matches organism names associated with reference genomes.[optional]")
target_args.add_argument("--tax-ids",   metavar = "TI-REGEX,", action = "store", help = "A valid regular expression that matches NCBI Taxonomy ID numbers for each genome to call against.[optional]")
target_args.add_argument("--gene-ids",  metavar = "GI-REGEX,", action = "store", help = "A valid regular expression that matches NCBI Gene ID numbers for each genome to call against.[optional]")
target_args.add_argument("--keep-all-sequences", action="store_true", default=False, help = "Do NOT discard any sequence in the database that is NOT a complete genome or complete plasmid sequence [optional]")

snp_caller_args = argparser.add_argument_group("SNP Caller Options")
snp_caller_args.add_argument('-s','--snp-caller', action="store", default = "samtools", choices = ["samtools"], help="Select the SNP Calling Program.[default:samtools]")
snp_caller_args.add_argument('--snp-caller-path', action="store", default = None, help = "Location of SNP Caller program binaries. [default to config setting]")
snp_caller_args.add_argument('--coverage', action="store_true", default=False, required=False, help = "Compute per-base coverage of alignment")

snp_anno_args = argparser.add_argument_group("Variant Annotation Options")
snp_anno_args.add_argument("--cache-dir", default = None, action = "store", type=str, help="The location to store raw and processed annotation source data. [default to config setting]")
snp_anno_args.add_argument('--snpeff-path', default = None, action = 'store', required = False, help = "Path to the snpEff.jar and .config files [default: default to config setting]")
snp_anno_args.add_argument('--blast-path', default = None, action = 'store', required = False, help = 'Path to the NCBI BLAST+ executables. [default: search system path]')

# snp_anno_args.add_argument('--blast-max', default = None, action = 'store', required = False, help = "Maximum score the number of Blast hits contribute to heuristic [default:20]")
# snp_anno_args.add_argument('--snp-max', default = None, action = 'store', required = False, help = "Maximum score the number of variants contribute to heuristic [default:20]")
snp_filt_args = argparser.add_argument_group("Variant Filtering Options")
for filter_type in EXPOSED_FILTERS:
	filter_type.customize_parser(snp_filt_args)

def main(args):

	configuration = get_config(args.config, alert = True)

	opts = {}
	opts['verbose'] = load_param(args.verbose, configuration['verbose'], False)
	opts['clean'] = load_param(args.clean, configuration["clean"], False)
	opts['cache_dir'] = load_param(args.cache_dir, configuration["cache_directory"], '.anno_cache')

	args.snp_caller_path = load_param(args.snp_caller_path, configuration["tool_paths"][args.snp_caller], "")

	args.snpeff_path = load_param(args.snpeff_path, configuration['tool_paths']['snpeff'], "")
	args.blast_path  = load_param(args.blast_path, configuration['tool_paths']['blast'], "")

	if not os.path.exists(args.sam_file): raise IOError("Input .sam File Not Found")
	start_clock = time()

	from pathovar.snp_caller.__main__ import call_snps

	call_result_files = call_snps(args, **opts)
	variant_file = call_result_files['variant_file']
	snp_called_time = time()
	if args.verbose: print('SNP Calling Done (%s sec)' % str(snp_called_time - start_clock))
	variant_locator_driver = None
	anno_vcf = None

	filter_args = utils.Namespace()
	filter_args.alt_depth =    load_param(args.alt_depth, configuration['filter_parameters']['alt_depth'], FilterByAltCallDepth.default_alt_depth)
	filter_args.min_depth =    load_param(args.min_depth, configuration['filter_parameters']['min_depth'], FilterByReadDepth.default_min_depth)
	filter_args.min_mq    =    load_param(args.min_mq, configuration['filter_parameters']['min_mq'],       FilterByMappingQuality.default_min_mq)
	filter_args.min_qual  =    load_param(args.min_qual, configuration['filter_parameters']['min_qual'],   FilterByCallQuality.default_min_qual)
	filter_args.ref_vcfs  =    load_param(args.ref_vcfs, configuration['filter_parameters']['ref_vcfs'],   FilterByComparisonVCF.default_ref_vcfs)
	filter_args.intersection = load_param(args.intersection, configuration['filter_parameters']['intersection'], FilterByComparisonVCF.default_intersection)

	from pathovar.snp_annotation.__main__ import run_annotation_report, find_variant_locations
	from pathovar.web.annotation_manager import EntrezAnnotationManager

	annotation_manager_driver = EntrezAnnotationManager(**opts)

	anno_vcf, variant_locator_driver = find_variant_locations(variant_file, annotation_manager_driver = annotation_manager_driver,
		**dict(filter_args = filter_args, **opts))

	annotation_report_driver = run_annotation_report(args, anno_vcf, variant_locator_driver, annotation_manager_driver, **opts)

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
