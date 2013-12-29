##@locate_variant
# Defines classes for mapping the location of SNPs into known genes 
# and annotated regions of genes.

## Standard Library Imports
from time import sleep
from collections import namedtuple

## Third Party Library Imports
import vcf 

## Internal Imports
from pathovar.web import entrez_eutils
from pathovar.snp_caller import snp_utils


##
#
class EntrezAnnotationMapper(object):

	##
	#
	def __init__(self, vcf_file, opts):
		self.vcf_file = vcf_file
		self.opts = opts
		self.verbose = opts["verbose"]
		self.entrez_handle = entrez_eutils.EntrezEUtilsDriver(self.opts)
		self.variants = [var for var in vcf.VCFReader(open(vcf_file, 'r'))]
		self.annotated_variants = []
		self.annotation_cache = dict()
	##
	#
	def annotate_snp(self, variant):
		ids = snp_utils.defline_parser(variant.CHROM)
		if len(ids.keys()) == 0:
			raise 
		key = ids.keys()[0]
		annotations = []
		if key is "gene_id":
			if 'gene_id-'+ids['gene_id'] not in self.annotation_cache:
				self.annotation_cache['gene_id-'+ids['gene_id']] = self.entrez_handle.find_nucleotide_annotations_by_gene_id(ids['gene_id'])
			annotations.extend(self.annotation_cache['gene_id-'+ids['gene_id']].locate_snp_site(variant.POS))
		annotation_records = map(self.fetch_gene_record, annotations)
		return annotation_records

	def find_nucleotide_annotations_by_gene_id(self, gene_id, mode = 'xml'):
		# Just in case it was passed as an int
		gene_id = str(gene_id)
		result = self.entrez_handle.find_nucleotides_by_gene_id(gene_id, form = 'genbank', mode = mode)
		if self.verbose:
			open(gene_id + '.genbank.' + mode, 'w').write(result)
		return GenBankFeatureFile(result, self.opts)

	def find_protein_annotations_by_gene_id(self, gene_id, mode = 'xml'):
		# Just in case it was passed as an int
		gene_id = str(gene_id)
		self.entrez_handle.find_protein_by_gene_id
		if self.verbose:
			open(gene_id + '.peptide_record.' + mode, 'w').write(get_protein_sequences_response.text)
		return GenBankFeatureFile( get_protein_sequences_response.text, self.opts)


	##
	#
	def fetch_gene_record(self, feature):
		if 'gene_id-'+feature.gid not in self.annotation_cache:
			sleep(1)
			self.annotation_cache['gene_id-'+feature.gid] = self.entrez_handle.find_protein_annotations_by_gene_id(feature.gid)
		anno_record = self.annotation_cache['gene_id-'+feature.gid]

		return feature

	def annotate_all_snps(self):
		for variant in self.variants:
			annotations = self.annotate_snp(variant)
			anno_variant = AnnotatedVariant(variant, annotations)
			if self.verbose:
				print(anno_variant)
			self.annotated_variants.append(anno_variant)

class AnnotatedVariant(vcf.model._Record):
	def __init__(self, _record, annotations = None):
		if annotations == None:
			annotations = []
		vcf.model._Record.__init__(self, _record.CHROM, _record.POS, _record.ID, _record.REF,
		 _record.ALT, _record.QUAL, _record.FILTER, _record.INFO, _record.FORMAT, 
		 _record._sample_indexes, _record.samples)
		self.annotations = annotations

	def __repr__(self):
		rep = "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s), ANNO=%(annotations)s" % self.__dict__
		return rep

	def __str__(self):
		rep = "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s), ANNO=%(annotations)s" % self.__dict__
		return rep

