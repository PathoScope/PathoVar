import os
import sys
import argparse
import subprocess
from time import time

import vcf

import pathovar
from pathovar import utils
from pathovar.utils.vcf_utils import EXPOSED_FILTERS
from pathovar.snp_annotation import snpeff_driver
from pathovar.snp_annotation.annotation_report import AnnotationReport
from pathovar.snp_annotation.locate_variant import VariantLocator
from pathovar.web.annotation_manager import EntrezAnnotationManager

argparser = argparse.ArgumentParser(prog = "snp_annotation")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("vcf_file", help = "The variant call file to annotate")
argparser.add_argument("--test", action = "store_true", required = False)
argparser.add_argument("--cache-dir", action = "store", type=str, default='.anno_cache', help="The location to store raw and processed annotation source data. [default='.anno_cache/']")
argparser.add_argument("--coverage", action = "store", type=str, default=False, help="The path to the coverage map .json file produced by calling snps with the --coverage setting." \
																							"If included, coverage results will be included in the output")
vcf_filter_args = argparser.add_argument_group("VCF Filters")
for filter_type in EXPOSED_FILTERS:
	filter_type.customize_parser(vcf_filter_args)
argparser.add_argument('--snpeff-path', default = '', action = 'store', required = False, help = "Path to the snpEff.jar and .config files [default: search system path]")
argparser.add_argument('--blast-path', default = '', action = 'store', required = False, help = 'Path to the NCBI BLAST executables. [default: search system path]')

def main(args):
	#if args.verbose: print(args)
	opts = dict()
	opts['verbose'] = args.verbose
	opts['cache_dir'] = args.cache_dir

	filter_args = utils.Namespace()
	filter_args.alt_depth = args.alt_depth
	filter_args.min_depth = args.min_depth
	filter_args.min_mq = args.min_mq
	filter_args.min_qual = args.min_qual
	filter_args.ref_vcfs = args.ref_vcfs
	filter_args.intersection = args.intersection

	opts['filter_args'] = filter_args
	coverage_data = None
	if args.coverage:
		from pathovar.snp_caller import compute_sam_coverage
		coverage_data = compute_sam_coverage.from_json(args.coverage)
	if not os.path.exists(args.vcf_file): raise IOError("Input .vcf File Not Found")
	timer = time()
	annotation_manager_driver = EntrezAnnotationManager(**opts)

	# Retrieve the annotated VCF and the VariantLocator instance for future reuse
	anno_vcf, variant_locator_driver = find_variant_locations(args.vcf_file, annotation_manager_driver, **opts)
	
	annotation_report_driver = run_annotation_report(args, anno_vcf, variant_locator_driver, 
		annotation_manager_driver, coverage_data = coverage_data, **opts)

	anno_json = annotation_report_driver.to_json_file()

	from pathovar.visualize.build_html_report import build_report
	build_report(anno_json)

	print("Annotation Complete (%r sec)" % (time() - timer))
	if(args.test):
		import IPython
		IPython.embed()

def find_variant_locations(vcf_file, annotation_manager_driver, **opts):
	variant_locator_driver = VariantLocator(vcf_file, annotation_manager = annotation_manager_driver, **opts)
	variant_locator_driver.annotate_all_snps()
	anno_vcf = variant_locator_driver.write_annotated_vcf()
	return anno_vcf, variant_locator_driver

def run_snpeff(args, anno_vcf, annotation_report_driver, **opts):
	try:
		if args.verbose: print("Running snpEff")
		if args.snpeff_path is None:
			raise Exception("snpEff path not set. Step Not Run.")
		eff_data = snpeff_driver.main(args.snpeff_path, anno_vcf, tempDir = None, 
			gidMap = annotation_report_driver.annotation_manager.genome_to_accesion_and_codon_table(), **opts)
		annotation_report_driver.consume_snpeff_results(eff_data)
	except snpeff_driver.snpEffException, e:
		print(e)

def run_annotation_report(args, anno_vcf, variant_locator_driver, annotation_manager_driver, coverage_data = None, **opts):
	
	annotation_report_driver = AnnotationReport(vcf_path=anno_vcf,variant_locator=variant_locator_driver, 
												annotation_manager = annotation_manager_driver, **opts)
	annotation_report_driver.merge_intergenic_record_chunks()
	#try:
	if args.coverage and coverage_data is not None:
		annotation_report_driver.compute_coverage_span(coverage_data)
	#except Exception, e:
	#	print("Error occurred computing coverage, %r" % e)
	ref_prot_fa = annotation_report_driver.generate_reference_protein_fasta_for_variants()

	# Load internal configuration file
	external_database_conf = pathovar.get_external_databases_config()
	enabled_databases = [database_name for database_name, database_conf in external_database_conf.items() if database_name != 'version' and  database_conf['enabled']]
	external_database_results = {}
	waiting_jobs = []
	
	# Path to configuration
	conf_path = pathovar.INSTALL_DIR
	try:
		# Opt-In Databases Job Queue
		if "comprehensive_antibiotic_resistance_database" in enabled_databases:
			from pathovar.snp_annotation.comprehensive_antibiotic_resistance_database_annotator import CARDProteinBlastAnnotator
			card_blast = CARDProteinBlastAnnotator(storage_path = os.path.join(conf_path, external_database_conf["comprehensive_antibiotic_resistance_database"]['storage_path']), bin_dir = args.blast_path)
			card_blast.query_with_proteins(ref_prot_fa)
			waiting_jobs.append(card_blast)


		if "drugbank" in enabled_databases:
			from pathovar.snp_annotation.drugbank_annotator import DrugBankProteinBlastAnnotator
			drugbank_blast = DrugBankProteinBlastAnnotator(storage_path = os.path.join(conf_path, external_database_conf['drugbank']['storage_path']), bin_dir = args.blast_path)
			drugbank_blast.query_with_proteins(ref_prot_fa)
			waiting_jobs.append(drugbank_blast)

		annotation_report_driver.get_entrez_gene_annotations()
		annotation_report_driver.get_entrez_biosystem_pathways()

		run_snpeff(args, anno_vcf, annotation_report_driver, **opts)

		# Block while annotations run
		for job in waiting_jobs:
			external_database_results[job.collection_name] = job.wait_for_results()

		# Consume the completed Blast searches 
		for external_database in external_database_results:
			for category in external_database_results[external_database]:
				annotation_report_driver.consume_blast_results(category, external_database_results[external_database][category])

		annotation_report_driver.normalize_all_entries()
		annotation_report_driver.score_all_entries()

	except ImportError, e:#Exception, e:
		print(e)

	return annotation_report_driver


if __name__ == '__main__':
	main(argparser.parse_args())
