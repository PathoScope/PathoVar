##@locate_variant
# Defines classes for mapping the location of SNPs into known genes 
# and annotated regions of genes.

## Standard Library Imports
import re
import os
import json
from time import sleep, time
from copy import copy
from collections import OrderedDict, namedtuple

## Third Party Library Imports
import vcf 
from vcf.parser import _Info
from bs4 import BeautifulSoup

## Internal Imports
from pathovar.web import entrez_eutils

from pathovar import utils
from pathovar.utils import vcf_utils

QUAL_FILTERS = [vcf_utils.FilterByAltCallDepth, vcf_utils.FilterByReadDepth]

def init_quality_filters(filter_args):
	filters = [filt(filter_args) for filt in QUAL_FILTERS]
	return filters

DEV_ANNOTATION_SCHEMA_STORE = dict()

##
#
class EntrezAnnotationMapper(object):
	##
	#
	def __init__(self, vcf_file, **opts):
		self.vcf_file = vcf_file
		self.opts = opts
		self.verbose = opts.get("verbose", False)
		self.cache_dir = opts.get("cache_dir", os.path.dirname(vcf_file) + os.sep + '.anno_cache')
		try:
			# If the cache directory is in the working directory, omit the leading os.sep
			# so that the system doesn't think it wants access to the root directory
			if self.cache_dir == os.sep + opts.get('cache_dir', '.anno_cache'):
				self.cache_dir = self.cache_dir[1:]
			os.makedirs(self.cache_dir)
		except OSError, e:
			pass
		self.no_cache = opts.get("no_cache", False)
		self.entrez_handle = entrez_eutils.EntrezEUtilsDriver(**opts)
		self.reader = vcf.VCFReader(open(vcf_file, 'r'))
		self.variants = vcf_utils.filter_vcf_in_memory(self.reader, init_quality_filters(self.opts['filter_args']), keep = True)
		self.annotated_variants = []
		self.annotation_cache = dict()

		self.reader.infos['GENE'] = _Info('GENE', 1, "String", "Gene containing this variant")

	##
	#
	def annotate_snp(self, variant):
		ids = utils.defline_parser(variant.CHROM)
		annotations = []
		opt_args = dict(verbose=self.verbose)
		if "gene_id" in ids:
			if 'gene_id-'+ids['gene_id'] not in self.annotation_cache:
				cache_file = self.cache_dir + os.sep + 'gene_id-'+ids['gene_id'] + '.annot.json'
				if self.verbose: print("Searching for cache file: %s" % cache_file)
				if os.path.exists(cache_file) and not self.no_cache:
					if self.verbose: print('Cache Hit')
					try:
						json_dict = json.load(open(cache_file))
						self.annotation_cache['gene_id-'+ids['gene_id']] = GenBankFeatureFile(json_dict, json=True, mol_type='nucl', **opt_args)
					except ValueError, e: 
						if self.verbose: print(e)
						self.annotation_cache['gene_id-'+ids['gene_id']] = self.find_annotations_by_gene_id(ids['gene_id'], xml=True, **opt_args)

				else:
					if self.verbose and not self.no_cache: print('Cache Miss')
					self.annotation_cache['gene_id-'+ids['gene_id']] = self.find_annotations_by_gene_id(ids['gene_id'], xml=True, **opt_args)
			if self.verbose: print("Localizing variant at positon %d" % variant.POS)
			annotations.extend(self.annotation_cache['gene_id-'+ids['gene_id']].locate_snp_site(variant.POS))
		return annotations

	##
	# Merge Nucleotide and Protein file annotations together for presentation to the 
	# user. Sequence Entries between Nucleotide and Protein files are keyed in the same
	# way, but the Protein file will contain different annotations.
	def find_annotations_by_gene_id(self, gene_id, **opts):
		opts['mol_type'] = 'nucl'
		nucleotide_annotations = self.find_nucleotide_annotations_by_gene_id(gene_id, **opts)
		opts['mol_type'] = 'prot'
		# The protein genpep file for the genome does not contain region annotations
		protein_annotations = self.find_protein_annotations_by_gene_id(nucleotide_annotations.features[0].gid, **opts)
		for gid, prot_entry in protein_annotations.entries.items():
			for annot_name, annot in prot_entry.annotations.items():
				nucleotide_annotations.entries[gid].annotations[annot_name].regions.extend(annot.regions)
		return nucleotide_annotations

	##
	#
	def find_nucleotide_annotations_by_gene_id(self, gene_id, **opts):
		# Just in case it was passed as an int
		gene_id = str(gene_id)
		result = self.entrez_handle.find_nucleotides_by_gene_id(gene_id, form = 'genbank', mode = 'xml')
		if self.verbose and not os.path.exists(gene_id + '.genbank.' + 'xml'):
			open(self.cache_dir + os.sep + gene_id + '.genbank.' + 'xml', 'w').write(result)
		return GenBankFeatureFile(result, **opts)
	##
	# 
	def find_protein_annotations_by_gene_id(self, gene_id, **opts):
		# Just in case it was passed as an int
		gene_id = str(gene_id)
		result = self.entrez_handle.find_protein_by_gene_id(gene_id, form = 'genbank', mode = 'xml')
		if self.verbose and not os.path.exists(gene_id + '.peptide_record.' + 'xml'):
			open(self.cache_dir + os.sep + gene_id + '.peptide_record.' + 'xml', 'w').write(result)
		return GenBankFeatureFile( result, **opts)

	##
	#
	def annotate_all_snps(self):
		for variant in self.variants:
			annotations = self.annotate_snp(variant)
			anno_variant = AnnotatedVariant(variant, annotations)
			self.annotated_variants.append(anno_variant)
		self.save_cache()

	def save_cache(self):
		if self.verbose: print("Saving Anntation Cache")
		for key,cached_annotator in self.annotation_cache.items():
			if cached_annotator.mol_type == 'nucl': open(self.cache_dir + os.sep + key + '.annot.json', 'w').write(json.dumps(cached_annotator.to_json_safe_dict()))

	def write_annotated_vcf(self):
		if self.verbose: print("Writing annotated .vcf file")
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

## GenBankFeatureFile
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the XML definition of a GenBank flat file.
class GenBankFeatureFile(object):
	def __init__(self, data, **opts):
		self.opts = opts
		self.verbose = opts['verbose']
		self.mol_type = opts['mol_type']
		timer = time()
		self.parser = None
		self.chromosome = None
		self.org_name = None
		self.entries = {}
		self.features = []
		if 'xml' in opts:
			self._parse_xml(data)
		elif 'json' in opts:
			self._from_json(data)

	def _parse_xml(self, xml):
		timer = time()
		self.parser = BeautifulSoup(xml)
		if self.verbose: print("XML Digested (%s sec)" % str(time() - timer))
		if self.verbose: print("Searching for Chromosome")
		self.org_name = self.parser.find("orgname_name").get_text().replace('\n','')
		self.chromosome = self.parser.find_all('iupacna')
		if self.chromosome:
			self.chromosome = ''.join(map(lambda x: x.get_text().strip(), self.chromosome))
			if self.verbose: print("Found, %d bp" % len(self.chromosome))

		if self.verbose: print("Gathering Entries and Features")
		self.entries = {ent.gid : ent for ent in map(lambda x: GenBankSeqEntry(x, self, xml = True), self.parser.find_all("seq-entry")) }
		self.features = map(lambda x: GenBankFeature(x, self), self.parser.find_all("seq-feat"))
		
		# Features are subsets of Entries. It would be a good idea to compress them to a single entity
		# later. Entries capture finer resolution details about a particular gene
		for feature in self.features:
			 if feature.gid in self.entries:
				 feature.title = self.entries[feature.gid].title
		
		# Keep only features that are not complete genomes
		self.features = [feature for feature in self.features if "complete genome" not in feature.title]
		
		# Each feature occurs multiple times in the file, redundant with its multiple regions. The complete
		# genomic span is the largest span. This works for single-span entities. 
		if self.verbose: print("Computing Genomic Coordinates")
		feature_dict = {}
		chromosome_len = len(self.chromosome)
		for feat in self.features:
			if feat.starts:
				if feat.gid not in feature_dict:
					feature_dict[feat.gid] = feat
				else:
					if feat.end - feat.start > feature_dict[feat.gid].end - feature_dict[feat.gid].start:
						feature_dict[feat.gid] = feat
					if feat.end > chromosome_len:
						print("%r exceeds chromosome size" % feat)
		self.features = feature_dict.values()
		
		# Using the genomic position data just computed, update coordinate information for the 
		# related entries
		for feat in self.features:
			entry = self.entries[feat.gid]
			entry.update_genome_position(feat)

	def _from_json(self, json_dict):
		if self.verbose: print("Starting to load from JSON")
		timer = time()
		self.org_name = json_dict['name']
		self.entries = {k:GenBankSeqEntry(v, self, **self.opts) for k,v in json_dict['entries'].items()}
		if self.verbose: print("Loading from JSON Complete (%rs)" % (time() - timer))


	def to_json_safe_dict(self):
		data_dict = {}
		data_dict["name"] = self.org_name
		data_dict['entries'] = {k: v.to_json_safe_dict() for k,v in self.entries.items() if "complete genome" not in v.title}
		return(data_dict)

	def locate_snp_site(self, snp_loc):
		last_entry_ind = 0
		entries = self.entries.values()
		for ind, entry in enumerate(entries[last_entry_ind:]):
			#if self.verbose: print(entry.start, entry.end)
			#if self.verbose: print("check %d >= %d and %d <= %d" % (snp_loc, entry.start, snp_loc, entry.end))
			if snp_loc >= entry.start and snp_loc <= entry.end:
				if self.verbose: print("SNP Location %s mapped within %r" % (snp_loc, entry))
				last_entry_ind += ind
				yield entry

	def __repr__(self):
		rep = "GenBankFile(" + self.mol_type + "|" + ', '.join(map(repr, self.entries)) + ")"
		return rep

## GenBankFeature
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the substructure of a GenBank flat file related to a single sequence feature
class GenBankFeature(object):
	def __init__(self, element, owner):
		self.element = element
		self.gid = element.find('seq-id_gi').get_text().strip()
		# Some features, especially in complex organisms, will have multi-part features. 
		# Capture all of that data
		self.starts = map(lambda x: int(x.get_text().strip()), self.element.find_all('seq-interval_from'))
		self.ends = map(lambda x: int(x.get_text().strip()), self.element.find_all('seq-interval_to'))

		# and set the extrema to the global start and stop
		self.start = None
		self.end   = None

		if(self.starts):
			self.start = min(self.starts)
			self.end = max(self.ends)
		else:
			#print("No Coordinates for %s" % self.gid)
			pass

		# Likely to be empty
		self.dbxref = element.find_all('dbtag_db')
		self.dbxref_type = map(lambda x: x.find("dbtag_db"), self.dbxref)
		self.dbxref_id = map(lambda x: x.find("object-id_id"), self.dbxref)

		self.title = None

	def __repr__(self):
		rep = "GenBankFeature(gi|%(gid)s, %(start)r, %(end)r, %(title)s)" % self.__dict__
		return rep

## GenBankSeqEntry
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the substructure of a GenBank flat file related to a single sequence entry
class GenBankSeqEntry(object):
	def __init__(self, data, owner, **opts):
		self.parser = None
		self.owner = owner

		self.gid = None
		self.accession = None
		
		# The coordinates in the seq-entry itself are relative to their own start and stop points. 
		# The genomic coordinates must be inferred from the genomic seq-feat tags
		self.starts = None
		self.ends   = None
		self.start  = None
		self.end	= None
		self.strand = None
		self.amino_acid_sequence = None
		# Prune later once starts and ends are set
		self.nucleotide_sequence = None 

		self.title = None
		self.annotations = {}
		self.comments = []
		
		self._id = None
		self._desc = None
		self._inst = None
		self._annot = None

		if 'xml' in opts:
			self._parse_xml(data)
		elif 'json' in opts:
			self._from_json(data)

	def _parse_xml(self, parser):
		self.parser = parser

		self.gid = parser.find('seq-id_gi').get_text().strip()
		self.accession = parser.find('textseq-id_accession').get_text().strip()
		
		self.strand = self.parser.find('na-strand')
		if self.strand:
			self.strand = self.strand.attrs['value']

		self.amino_acid_sequence = map(lambda x: x.get_text(), parser.find_all('iupacaa'))
		# Prune later once starts and ends are set
		self.nucleotide_sequence = self.owner.chromosome

		self.title = parser.find('seqdesc_title').get_text().strip()
		self.annotations = {ann.name: ann for ann in map(lambda d: GenBankAnnotation(d, xml=True), parser.find_all("bioseq_annot"))}
		self.comments = map(lambda x: x.get_text().strip(), parser.find_all('seq-feat_comment'))
		
		#self._id = parser.find_all('bioseq_id')
		#self._desc = parser.find_all("bioseq_descr")
		#self._inst = parser.find_all("bioseq_inst")
		#self._annot = parser.find_all("bioseq_annot")


	def _from_json(self, json_dict):
		self.gid = json_dict['gid']
		self.accession = json_dict['accession']
		
		self.strand = json_dict['strand']
		self.amino_acid_sequence = json_dict['amino_acid_sequence']
		self.nucleotide_sequence = json_dict['nucleotide_sequence']
		
		self.title = json_dict['title']
		self.comments = json_dict['comments']
		
		self.start =  json_dict["start"]
		self.end = json_dict["end"]
		self.starts = json_dict["starts"]
		self.ends =  json_dict["ends"]
		
		self.annotations = {k:GenBankAnnotation(v, json=True) for k,v in json_dict['annotations'].items()}

	def update_genome_position(self, feature):
		self.starts = feature.starts
		self.ends   = feature.ends  
		self.start  = feature.start 
		self.end    = feature.end
		if self.owner.mol_type == 'nucl':
			self.nucleotide_sequence = self.owner.chromosome[self.start:self.end]

	def __repr__(self):
		rep = "GenBankSeqEntry(gi=%(gid)s, Title=%(title)s, ACC=%(accession)s)" 
		return rep % self.__dict__

	def to_info_field(self):
		rep = "(gi:%(gid)s|title:%(title)s|acc:%(accession)s|strand:%(strand)s" % self.__dict__
		annos = map(lambda x: x.to_info_field(), self.annotations.values())
		rep += '||Annotations|' + '|'.join(annos)
		rep += '||Comments|' + '|'.join(self.comments)
		if re.search(r'resistance', rep):
			rep += '||RESISTANCE'
		rep += ')'
		rep = re.sub(r'\s', '_', rep)
		return rep

	def to_json_safe_dict(self):
		data_dict = copy(self.__dict__)
		# Remove fields that don't serialize.
		unsafe_fields = ["parser", "_id","_desc","_inst","_annot","owner"]
		for f in unsafe_fields:
			try:
				data_dict.pop(f)
			except KeyError:
				pass

		# Translate nested objects to dictionaries
		data_dict["annotations"] = {k:v.__dict__ for k,v in data_dict["annotations"].items()}

		return data_dict

class AnnotationExtension(OrderedDict):
	def __init__(self,*args, **kwargs):
		OrderedDict.__init__(self, *args, **kwargs)

	def __repr__(self):
		rep = ''
		ext_type = self.get("ext_type", 'anno_ext')
		rep += ext_type + '('
		vals = []
		for key in self:
			if key == "ext_type": continue
			vals.append(key + ': ' + str(self[key]))
		rep += ', '.join(vals) + ')'
		return rep

	def __str__(self):
		return repr(self)

class GenBankAnnotation(object):
	def __init__(self, data, **opts):
		self.start = []
		self.end = []
		self.name = None
		self.regions = []

		if 'xml' in opts:
			self._parse_xml(data)
		elif 'json' in opts:
			self._from_json(data)

	def _parse_xml(self, data):
		self.start = map(lambda x: int(x.get_text().replace(' ','')), data.find_all("seq-interval_from"))
		self.end = map(lambda x: int(x.get_text().replace(' ','')), data.find_all("seq-interval_to"))
		
		self.regions = [] #list(map(lambda x: x.get_text().strip(), data.find_all("seqfeatdata_region")))
		self.name = ', '.join(map(lambda x: x.get_text().strip(), data.find_all("prot-ref_name_e")))
		raw_extensions = data.find_all("seq-feat_ext")
		for raw_ext in raw_extensions:
			obj_type = str(raw_ext.find("object-id_str").get_text().strip())
			obj_data = map(lambda x: x.get_text().strip(), raw_ext.find_all("user-object_data")[0].find_all("user-field_data"))
			ext = AnnotationExtension()
			DEV_ANNOTATION_SCHEMA_STORE[obj_type] = raw_ext
			ext = self.process_extension(raw_ext)
			ext['ext_type'] = obj_type
			self.regions.append(ext)

	def _from_json(self, json_dict):
		self.start = json_dict['start']
		self.end = json_dict['end']
		self.name = json_dict['name']
		self.regions = json_dict['regions']

	def __repr__(self):
		rep = "GenBankAnnotation(Starts=%(start)s, Ends=%(end)s, name:%(name)s, regions:%(regions)s)" % self.__dict__
		return rep

	def process_extension(self, raw_ext):
		labels = map(lambda x: x.get_text().strip(), raw_ext.find_all("user-field_label"))
		data = raw_ext.find_all("user-field_data")
		values = []
		for datum in data:
			text = datum.get_text().strip()
			if len(text) == 0:
				text = str(datum.attrs)
				if text == '{}':
					text = ''
			values.append(text)
		ext = AnnotationExtension()
		for i in xrange(len(labels)):
			ext[labels[i]] = values[i]

		return ext

	def to_info_field(self):
		rep = '(starts:%(start)s|ends:%(end)s|' % self.__dict__
		for region in self.regions:
			buff = []
			region = copy(region)
			ext_type = region.pop('ext_type')
			pos_from = region.pop('from', None)
			pos_to   = region.pop('to', None)
			for key, value in region.items():
				buff.append(key + ":" + str(value))
			# prepend indices
			if pos_from:
				buff = ['from:' + str(pos_from).strip(), 'to:' + str(pos_to).strip()] + buff
			buff = '|'.join(buff)
			buff = ext_type + '(' + buff
			rep += buff
			rep += ')'
		rep += ')'
		return rep

class UnknownAnnotationException(Exception):
	pass
