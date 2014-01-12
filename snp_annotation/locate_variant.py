##@locate_variant
# Defines classes for mapping the location of SNPs into known genes 
# and annotated regions of genes.

## Standard Library Imports
import json
import re
from time import sleep
from collections import OrderedDict

## Third Party Library Imports
import vcf 
from vcf.parser import _Info
from bs4 import BeautifulSoup

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
		self.reader = vcf.VCFReader(open(vcf_file, 'r'))
		self.variants = [var for var in self.reader]
		self.annotated_variants = []
		self.annotation_cache = dict()

		self.reader.infos['GENE'] = _Info('GENE', 1, "String", "Gene containing this variant")

	##
	#
	def annotate_snp(self, variant):
		ids = snp_utils.defline_parser(variant.CHROM)
		annotations = []
		opt_args = dict(verbose=self.verbose)
		if "gene_id" in ids:
			if 'gene_id-'+ids['gene_id'] not in self.annotation_cache:
				self.annotation_cache['gene_id-'+ids['gene_id']] = self.find_annotations_by_gene_id(ids['gene_id'], opt_args)
			annotations.extend(self.annotation_cache['gene_id-'+ids['gene_id']].locate_snp_site(variant.POS))
		return annotations


	##
	# Merge Nucleotide and Protein file annotations together for presentation to the 
	# user. Sequence Entries between Nucleotide and Protein files are keyed in the same
	# way, but the Protein file will contain different annotations.
	def find_annotations_by_gene_id(self, gene_id, opts):
		opts['mol_type'] = 'nucl'
		nucleotide_annotations = self.find_nucleotide_annotations_by_gene_id(gene_id, opts)
		opts['mol_type'] = 'prot'
		# The protein genpep file for the genome does not contain region annotations
		protein_annotations = self.find_protein_annotations_by_gene_id(nucleotide_annotations.features[0].gid, opts)
		for gid, prot_entry in protein_annotations.entries.items():
			for annot_name, annot in prot_entry.annotations.items():
				nucleotide_annotations.entries[gid].annotations[annot_name].regions.extend(annot.regions)
		return nucleotide_annotations
	##
	#
	def find_nucleotide_annotations_by_gene_id(self, gene_id, opts):
		# Just in case it was passed as an int
		gene_id = str(gene_id)
		result = self.entrez_handle.find_nucleotides_by_gene_id(gene_id, form = 'genbank', mode = 'xml')
		if self.verbose:
			open(gene_id + '.genbank.' + 'xml', 'w').write(result)
		return GenBankFeatureFile(result, opts)
	##
	# 
	def find_protein_annotations_by_gene_id(self, gene_id, opts):
		# Just in case it was passed as an int
		gene_id = str(gene_id)
		result = self.entrez_handle.find_protein_by_gene_id(gene_id, form = 'genbank', mode = 'xml')
		if self.verbose:
			open(gene_id + '.peptide_record.' + 'xml', 'w').write(result)
		return GenBankFeatureFile( result, opts)

	##
	#
	def annotate_all_snps(self):
		for variant in self.variants:
			annotations = self.annotate_snp(variant)
			anno_variant = AnnotatedVariant(variant, annotations)
			if self.verbose:
				print(anno_variant)
			self.annotated_variants.append(anno_variant)


	def write_annotated_vcf(self):
		output_file = self.vcf_file[:-4] + '.anno.vcf'
		writer = vcf.Writer(open(output_file, 'w'), self.reader)
		for variant in self.annotated_variants:
			writer.write_record(variant)
		writer.close()
		return output_file

##
#
class AnnotatedVariant(vcf.model._Record):
	def __init__(self, _record, annotations = None):
		if annotations == None:
			annotations = []
		vcf.model._Record.__init__(self, _record.CHROM, _record.POS, _record.ID, _record.REF,
		 _record.ALT, _record.QUAL, _record.FILTER, _record.INFO, _record.FORMAT, 
		 _record._sample_indexes, _record.samples)
		self.INFO['GENE'] = map(lambda x: x.to_info_field(), annotations)
		self.annotations = annotations

	def __repr__(self):
		rep = "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s), ANNO=%(annotations)s" % self.__dict__
		return rep

	def __str__(self):
		rep = "Record(CHROM=%(CHROM)s, POS=%(POS)s, REF=%(REF)s, ALT=%(ALT)s), ANNO=%(annotations)s" % self.__dict__
		return rep

# VariantTuple = namedtuple("VariantTuple",["alt_allele", "position"])
# class MutantSequenceFactory(object):
# 	def __init__(self, gid, name, reference_sequence, variant_list = None, opts = None):
# 		if variant_list is None:
# 			variant_list = []
# 		if opts is None:
# 			opts = dict(verbose = False)
# 		self.gid = gid
# 		self.name = name
# 		self.reference_sequence = reference_sequence
# 		self.variants = variant_list
# 		self.opts = opts
# 		self.mutant_sequences = []

# 	def compute_mutants(self):
# 		for variant in self.variants:
# 			for alt in variant.ALT:
# 				variant_tuple = VariantTuple(variant.alt, variant.POS, [])


# class MutantSequence(list):
# 	## 
# 	# variant_positions is a list of 3-tuples, (ALT_ALLELE, POSITION, MUTATION_TYPE_LIST)
# 	def __init__(self, reference_sequence, variant_positions, annotations):
# 		list.__init__(self, reference_sequence)
# 		self.variant_positions = variant_positions
# 		self.annotations = annotations
# 		for variant in self.variant_positions:
# 			self[variant[1]] = variant[0]
# 			for annotation in self.annotations:
# 				if annotation.start <= variant[1] <= annotation.end:
# 					variant[2].append(annotation)
# 					for region in annotation.regions:
# 						if region.start <= variant[1] <= region.end:
# 							variant[2].append(region)




## GenBankFeatureFile
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the XML definition of a GenBank flat file.
class GenBankFeatureFile(object):
	def __init__(self, xml, opts):
		self.parser = BeautifulSoup(xml)
		self.opts = opts
		self.verbose = opts['verbose']
		self.mol_type = opts['mol_type']
		self.chromosome = self.parser.find('iupacna')
		if self.chromosome:
			self.chromosome = self.chromosome.get_text()
		self.entries = {ent.gid : ent for ent in map(lambda x: GenBankSeqEntry(x, self), self.parser.find_all("seq-entry")) }
		self.features = map(GenBankFeature, self.parser.find_all("seq-feat"))
		
		# Features are subsets of Entries. It would be a good idea to compress them to a single entity
		# later. Entries capture finer resolution details about a particular gene
		for feature in self.features:
			 if feature.gid in self.entries:
				 feature.title = self.entries[feature.gid].title
		
		# Keep only features that are not complete genomes
		self.features = [feature for feature in self.features if "complete genome" not in feature.title]
		
		# Each feature occurs multiple times in the file, redundant with its multiple regions. The complete
		# genomic span is the largest span. This works for single-span entities. 
		feature_dict = {}
		for feat in self.features:
			if feat.gid not in feature_dict:
				feature_dict[feat.gid] = feat
			else:
				if feat.end - feat.start > feature_dict[feat.gid].end - feature_dict[feat.gid].start:
					feature_dict[feat.gid] = feat
		self.features = feature_dict.values()
		
		# Using the genomic position data just computed, update coordinate information for the 
		# related entries
		for feat in self.features:
			entry = self.entries[feat.gid]
			entry.update_genome_position(feat)


	def locate_snp_site(self, snp_loc):
		contains = []
		last_end = -1
		for entry in self.entries.values():
			if snp_loc >= entry.start and snp_loc <= entry.end:
				contains.append(entry)
				if self.verbose: print("SNP Location %s mapped within %r" % (snp_loc, entry))

		return contains

	def __repr__(self):
		rep = "(" + self.mol_type + "|" + ', '.join(map(repr, self.entries)) + ")"
		return rep

## GenBankFeature
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the substructure of a GenBank flat file related to a single sequence feature
class GenBankFeature(object):
	def __init__(self, element):
		self.element = element
		self.gid = element.find('seq-id_gi').get_text()
		# Some features, especially in complex organisms, will have multi-part features. 
		# Capture all of that data
		self.starts = map(lambda x: int(x.get_text()), self.element.find_all('seq-interval_from'))
		self.ends = map(lambda x: int(x.get_text()), self.element.find_all('seq-interval_to'))

		# and set the extrema to the global start and stop
		self.start = min(self.starts)
		self.end = max(self.ends)

		# Likely to be empty
		self.dbxref = element.find_all('dbtag_db')
		self.dbxref_type = map(lambda x: x.find("dbtag_db"), self.dbxref)
		self.dbxref_id = map(lambda x: x.find("object-id_id"), self.dbxref)

		self.title = None

	def __repr__(self):
		rep = "(gi:%(gid)s, %(start)r, %(end)r, %(title)s)" % self.__dict__
		return rep

## GenBankSeqEntry
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the substructure of a GenBank flat file related to a single sequence entry
class GenBankSeqEntry(object):
	def __init__(self, element, owner):
		self.element = element
		self.owner = owner

		self.gid = element.find('seq-id_gi').get_text()
		self.accession = element.find('textseq-id_accession').get_text()
		
		# The coordinates in the seq-entry itself are relative to their own start and stop points. 
		# The genomic coordinates must be inferred from the genomic seq-feat tags
		self.starts = None
		self.ends   = None
		self.start  = None
		self.end	= None
		
		self.amino_acid_sequence = element.find('iupacaa')
		# Prune later once starts and ends are set
		self.nucleotide_sequence = owner.chromosome

		self.title = element.find('seqdesc_title').get_text()
		self.annotations = {ann.name: ann for ann in map(GenBankAnnotation, element.find_all("bioseq_annot"))}
		
		self._id = element.find_all('bioseq_id')
		self._desc = element.find_all("bioseq_descr")
		self._inst = element.find_all("bioseq_inst")
		self._annot = element.find_all("bioseq_annot")

	def update_genome_position(self, feature):
		self.starts = feature.starts
		self.ends   = feature.ends  
		self.start  = feature.start 
		self.end	= feature.end   

		self.nucleotide_sequence = self.nucleotide_sequence[self.start:self.end]

	def __repr__(self):
		rep = "(gi:%(gid)s, Title:%(title)s, ACC:%(accession)s, ANNO:%(annotations)s)" % self.__dict__
		return rep

	def to_json(self):
		ref_dict = OrderedDict()
		ref_dict["gid"] = self.gid
		ref_dict["accession"] = self.accession
		ref_dict["starts"] = self.starts
		ref_dict["ends"] = self.ends
		ref_dict["title"] = self.title
		ref_dict["annotations"] = {k:v.__dict__ for k,v in self.annotations.items()}
		return json.dumps(ref_dict)

	def to_info_field(self):
		rep = "(%(gid)s|%(title)s|%(accession)s)" % self.__dict__
		rep = re.sub(r'[^a-zA-Z0-9_\(\)\|]', '_', rep)
		return rep

class GenBankAnnotation(object):
	def __init__(self, element):
		self.start = map(lambda x: int(x.get_text()), element.find_all("seq-interval_from"))
		self.end = map(lambda x: int(x.get_text()), element.find_all("seq-interval_to"))
		
		self.regions = [] #list(map(lambda x: x.get_text(), element.find_all("seqfeatdata_region")))
		self.name = ', '.join(map(lambda x: x.get_text(), element.find_all("prot-ref_name_e")))
		raw_extensions = element.find_all("seq-feat_ext")
		for raw_ext in raw_extensions:
			obj_type = str(raw_ext.find("object-id_str").get_text())
			obj_data = map(lambda x: x.get_text(), raw_ext.find_all("user-object_data")[0].find_all("user-field_data"))
			ext = OrderedDict()
			if obj_type == "cddScoreData":
				ext['start'] = int(obj_data[0])
				ext['end'] = int(obj_data[1])
				ext['definition'] = str(obj_data[2]).replace('\n','')
				ext['name'] = str(obj_data[3]).replace('\n','')
				ext['e_value'] = float(obj_data[5])
				ext['type'] = 'GenBankCDD'
				self.regions.append(ext)
			else:
				raise UnknownAnnotationException("Unknown Extension: %s" % obj_type)

	def __repr__(self):
		rep = "(starts:%(start)s, ends:%(end)s, name:%(name)s, regions:%(regions)s)" % self.__dict__
		return rep

class UnknownAnnotationException(Exception):
	pass