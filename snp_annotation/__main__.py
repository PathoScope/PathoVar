import os
import sys
import argparse
import subprocess
from time import time

import vcf

import pathovar
from pathovar import utils
from pathovar.snp_annotation import locate_variant, annotation_report
from pathovar.web import annotation_manager

argparser = argparse.ArgumentParser(prog = "snp_annotation")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("vcf_file", help = "The variant call file to annotate")
argparser.add_argument("--test", action = "store_true", required = False)
argparser.add_argument("--cache-dir", action = "store", type=str, default='.anno_cache', help="The location to store raw and processed annotation source data. [default='.anno_cache/']")
argparser.add_argument('--min-depth', type=int, default=5, help="The minimum number of reads that must map to a location to trust a given variant call [default:5]")
argparser.add_argument('--alt-depth', type=float, default=0.4, help="The MAF threshold, under which variants are ignored [default:0.4]")

def main():
	args = argparser.parse_args()
	#if args.verbose: print(args)
	opts = dict()
	opts['verbose'] = args.verbose
	opts['cache_dir'] = args.cache_dir

	filter_args = utils.Namespace()
	filter_args.alt_depth = args.alt_depth
	filter_args.min_depth = args.min_depth

	opts['filter_args'] = filter_args

	if not os.path.exists(args.vcf_file): raise IOError("Input .vcf File Not Found")
	timer = time()
	annotation_manager_driver = annotation_manager.EntrezAnnotationManager(**opts)
	variant_locator = locate_variant.VariantLocator(args.vcf_file, annotation_manager = annotation_manager_driver, **opts)
	variant_locator.annotate_all_snps()

	anno_vcf = variant_locator.write_annotated_vcf()

	annotation_report_driver = annotation_report.AnnotationReport(anno_vcf, annotation_manager_driver, **opts)
	ref_prot_fa = annotation_report_driver.generate_reference_protein_fasta_for_variants()

	# Load internal configuration file
	external_database_conf = pathovar.get_external_databases_config()
	enabled_databases = [database_name for database_name, database_conf in external_database_conf.items() if database_name != 'version' and  database_conf['enabled']]
	external_database_results = {}
	waiting_jobs = []
	
	# Opt-In Databases Job Queue
	if "comprehensive_antibiotic_resistance_database" in enabled_databases:
		from pathovar.snp_annotation.comprehensive_antibiotic_resistance_database_annotator import CARDProteinBlastAnnotator
		card_blast = CARDProteinBlastAnnotator()
		card_blast.query_with_proteins(ref_prot_fa)
		waiting_jobs.append(card_blast)


	if "drugbank" in enabled_databases:
		from pathovar.snp_annotation.drugbank_annotator import DrugBankProteinBlastAnnotator
		drugbank_blast = DrugBankProteinBlastAnnotator()
		drugbank_blast.query_with_proteins(ref_prot_fa)
		waiting_jobs.append(drugbank_blast)

	annotation_report_driver.get_entrez_gene_annotations()
	annotation_report_driver.get_entrez_biosystem_pathways()
	annotation_report_driver.merge_intergenic_record_chunks()

	# Block while annotations run
	for job in waiting_jobs:
		external_database_results[job.collection_name] = job.wait_for_results()

	# Consume the completed Blast searches 
	for external_database in external_database_results:
		for category in external_database_results[external_database]:
			annotation_report_driver.consume_blast_results(category, external_database_results[external_database][category])

	annotation_report_driver.to_json_file()

	print("Annotation Complete (%r sec)" % (time() - timer))
	if(args.test):
		import IPython
		IPython.embed()


if __name__ == '__main__':
	main()