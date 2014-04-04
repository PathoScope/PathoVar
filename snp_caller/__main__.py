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
snp_caller_args.add_argument('--snp-caller-path', action="store", default = "", help = "Location of SNP Caller program binaries. Default will search for them on the system path")
snp_caller_args.add_argument('-c','--coverage', action="store_true", default=False, required=False, help = "Compute per-base coverage of alignment")

def call_snps(args, **opts):
    snp_caller_driver = None
    if args.snp_caller == "samtools":
        from pathovar.snp_caller import samtools_snp_caller
        snp_caller_driver = samtools_snp_caller.SamtoolsSNPCaller(bin_dir = args.snp_caller_path, **opts)
    variant_file = snp_caller_driver.call_snps(args.sam_file, source = args.reference_genomes, 
        org_names_reg = args.org_names, tax_ids_reg = args.tax_ids, gene_ids_reg = args.gene_ids, 
        keep_all = args.keep_all_sequences)
    consensus_sequences = variant_file + ".cns.fq"
    result_files = {'variant_file': variant_file, 'consensus': consensus_sequences}
    if args.coverage:
        from pathovar.snp_caller import compute_sam_coverage
        coverage_json = compute_sam_coverage.main(sam_parser = snp_caller_driver.sam_parser)
        result_files['coverage'] = coverage_json
    return result_files

def main(args):
    opts = {}
    opts['verbose'] = args.verbose
    opts['clean'] = args.clean

    if not os.path.exists(args.sam_file): raise IOError("Input .sam File Not Found")
    start_clock = time()
    variant_file = call_snps(args, **opts)
    snp_called_time = time()
    if args.verbose: print('SNP Calling Done (%s sec)' % str(snp_called_time - start_clock))

if __name__ == '__main__':
    args = argparser.parse_args()
    main(args)