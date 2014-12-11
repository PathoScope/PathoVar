import os
import argparse
from time import time

from pathovar import utils
from pathovar.utils.config import load_param, get_config

from pathovar.utils.vcf_utils import EXPOSED_FILTERS
from pathovar.utils.vcf_utils import FilterByAltCallDepth
from pathovar.utils.vcf_utils import FilterByReadDepth
from pathovar.utils.vcf_utils import FilterByMappingQuality
from pathovar.utils.vcf_utils import FilterByCallQuality
from pathovar.utils.vcf_utils import FilterByComparisonVCF

from pathovar.snp_annotation import snpeff_driver
from pathovar.snp_annotation.blast_driver import BlastDriverException
from pathovar.snp_annotation.annotation_report import AnnotationReport
from pathovar.snp_annotation.locate_variant import VariantLocator
from pathovar.web.annotation_manager import EntrezAnnotationManager

argparser = argparse.ArgumentParser(prog = "snp_annotation")
argparser.add_argument("-v", "--verbose", action = "store_true", required = False)
argparser.add_argument("vcf_file", help = "The variant call file to annotate")
argparser.add_argument("--test", action = "store_true", required = False)
argparser.add_argument("-c", "--config", action = "store", default = None)
argparser.add_argument("--cache-dir", action = "store", type=str, default=None, help="The location to store raw and processed annotation source data. [default=]")
argparser.add_argument("--coverage", action = "store", type=str, default=False, help="The path to the coverage map .json file produced by calling snps with the --coverage setting." \
																							"If included, coverage results will be included in the output")
argparser.add_argument("--clean", action = "store_true", required = False, help = "Remove intermediary files after finished")
vcf_filter_args = argparser.add_argument_group("VCF Filters")
for filter_type in EXPOSED_FILTERS:
	filter_type.customize_parser(vcf_filter_args)
argparser.add_argument('--snpeff-path', default = None, action = 'store', required = False, help = "Path to the snpEff.jar and .config files [default: search system path]")
argparser.add_argument('--blast-path', default = None, action = 'store', required = False, help = 'Path to the NCBI BLAST executables. [default: search system path]')

CatchError = Exception
CatchError = ArithmeticError

def main(args):

	configuration = get_config(args.config, alert = True)

	#if args.verbose: print(args)
	opts = dict()
	opts['verbose'] = load_param(args.verbose, configuration['verbose'], False)
	opts['clean'] = load_param(args.clean, configuration["clean"], False)
	opts['cache_dir'] = load_param(args.cache_dir, configuration["cache_directory"], '.anno_cache')

	filter_args = utils.Namespace()
	filter_args.alt_depth =    load_param(args.alt_depth, configuration['filter_parameters']['alt_depth'], FilterByAltCallDepth.default_alt_depth)
	filter_args.min_depth =    load_param(args.min_depth, configuration['filter_parameters']['min_depth'], FilterByReadDepth.default_min_depth)
	filter_args.min_mq    =    load_param(args.min_mq, configuration['filter_parameters']['min_mq'],       FilterByMappingQuality.default_min_mq)
	filter_args.min_qual  =    load_param(args.min_qual, configuration['filter_parameters']['min_qual'],   FilterByCallQuality.default_min_qual)
	filter_args.ref_vcfs  =    load_param(args.ref_vcfs, configuration['filter_parameters']['ref_vcfs'],   FilterByComparisonVCF.default_ref_vcfs)
	filter_args.intersection = load_param(args.intersection, configuration['filter_parameters']['intersection'], FilterByComparisonVCF.default_intersection)

	args.snpeff_path = load_param(args.snpeff_path, configuration['tool_paths']['snpeff'], "")
	args.blast_path  = load_param(args.blast_path, configuration['tool_paths']['blast'], "")


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
		print("An error occurred while running SnpEff: %s" % e)

def run_annotation_report(args, anno_vcf, variant_locator_driver, annotation_manager_driver, coverage_data = None, **opts):
	annotation_report_driver = AnnotationReport(vcf_path=anno_vcf,variant_locator=variant_locator_driver, 
												annotation_manager = annotation_manager_driver, **opts)
	annotation_report_driver.merge_intergenic_record_chunks()
	try:
		if args.coverage and coverage_data is not None:
			annotation_report_driver.compute_coverage_span(coverage_data)
	except Exception, e:
		print("Error occurred while computing coverage, %r" % e)
	ref_prot_fa = annotation_report_driver.generate_reference_protein_fasta_for_variants()

	# Load internal configuration file
	conf = get_config(args.config)
	external_database_conf = conf['external_databases']
	enabled_databases = [database_name for database_name, database_conf in external_database_conf.items() if database_conf['enabled']]
	external_database_results = {}
	waiting_jobs = []

	try:
		# Opt-In Databases Job Queue
		if "comprehensive_antibiotic_resistance_database" in enabled_databases:
			from pathovar.snp_annotation.comprehensive_antibiotic_resistance_database_annotator import CARDProteinBlastAnnotator
			storage_path = (conf['database_storage_directory'] + os.sep + external_database_conf["comprehensive_antibiotic_resistance_database"]['storage_path'])
			card_blast = CARDProteinBlastAnnotator(storage_path = storage_path,
													bin_dir = args.blast_path, 
													clean = args.clean)
			card_blast.query_with_proteins(ref_prot_fa)
			waiting_jobs.append(card_blast)

		if "drugbank" in enabled_databases:
			from pathovar.snp_annotation.drugbank_annotator import DrugBankProteinBlastAnnotator
			storage_path = (conf['database_storage_directory']  + os.sep +  external_database_conf["drugbank"]['storage_path'])
			drugbank_blast = DrugBankProteinBlastAnnotator(storage_path = storage_path, 
															bin_dir = args.blast_path, 
															clean = args.clean)
			drugbank_blast.query_with_proteins(ref_prot_fa)
			waiting_jobs.append(drugbank_blast)			

		annotation_report_driver.get_entrez_gene_annotations()
		annotation_report_driver.get_entrez_biosystem_pathways()

		run_snpeff(args, anno_vcf, annotation_report_driver, **opts)
		# Block while annotations run
		for job in waiting_jobs:
			try:
				external_database_results[job.collection_name] = job.wait_for_results()
			except BlastDriverException, e:
				print("Waiting error", e)

		# Consume the completed Blast searches 
		for external_database in external_database_results:
			for category in external_database_results[external_database]:
				annotation_report_driver.consume_blast_results(category, external_database_results[external_database][category])

		# Do any cleaning if needed
		[job.clean() for job in waiting_jobs]

		annotation_report_driver.normalize_all_entries()
		annotation_report_driver.score_all_entries()

	except CatchError, e:
		print("An error occured while building report, %s" % e)

	return annotation_report_driver


if __name__ == '__main__':
	main(argparser.parse_args())
