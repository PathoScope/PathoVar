##@package entrez_eutils
# Defines logic for interacting with Entrez's EUtils web services. Used to fetch genome 
# and sequence annotations. 
#
# TODO:
# - Migrate All components of the GenBankFeatureFile
# system to crossreference_snp_location.py in snp_annotation/
# - Include sequence data in GenBankFeature for first approximation
# of mutation prediction

# System Dependencies
import sys
import re
import os
from collections import namedtuple

# External Dependencies
import requests
from bs4 import BeautifulSoup

## URL CONSTANTS

taxonomy_summary_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={tid}'
taxonomy_detail_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={tid}'
# Organism name to Taxonomy ID: 
org_name_to_taxonomy_id = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={org_name}[organism]"

# Genome Link by Org Name
genome_by_org_name = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term={org_name}[organism]'

# Set remote environment for translating from genome to nucleotides
link_from_genome_to_nuccore = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id={id}&cmd=neighbor_history'

# Retrieve all matching nucleotide sequences
get_nucleotides_from_link = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&' \
                            'query_key={query_key}&WebEnv={web_env}&rettype={form}&retmode={mode}'
# Retrieve a particular sequence record by gene id
get_nucleotides_by_gene_id = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={gene_id}&rettype={form}&retmode={mode}'

# Retrieve a particular sequence record by gene id
get_protein_by_gene_id = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={gene_id}&rettype={form}&retmode={mode}'

get_gene_by_gene_id = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={gene_id}&rettype={form}"

## EntrezEUtilsDriver
# Defines all of the logic for interacting with Entrez's EUtils web service, handling errors, 
# and performing multi-request actions.
#
class EntrezEUtilsDriver(object):
    def __init__(self, opts):
        self.opts = opts
        self.verbose = self.opts['verbose']

    def find_genome_by_org_name(self, org_name, form = 'fasta', mode = 'text'):
        # Find the genome id
        if self.verbose: print("Fetching Genome ID from Entrez")
        genome_db_response = requests.get(genome_by_org_name.format(**dict(org_name = org_name)))
        genome_db_response.raise_for_status()
        
        # Parse the response
        if self.verbose: print("Response recieved...")
        genome_db_response_xml = BeautifulSoup(genome_db_response.text)
        genome_id = genome_db_response_xml.find_all('id')
        if(len(genome_id) == 0):
            if self.verbose: print("No Genome ID found")
            output_message = genome_db_response_xml.find_all('outputmessage')[0]
            if output_message.get_text() == u"No items found.":
                raise OrganismNotFoundException()
        # If there is ambiguity over which genome id to choose, just take the first
        genome_id = genome_id[0].get_text()
        if self.verbose: print("Genome ID: %s" % genome_id)

        # Set up remote environment to cross-link from genome id to nuccore id
        if self.verbose: print("Fetching Link to Genome from Entrez")
        genome_to_nuccore_response = requests.get(link_from_genome_to_nuccore.format(**dict(id = genome_id)))
        genome_to_nuccore_response.raise_for_status()

        if self.verbose: print("Response recieved...")
        query_key = None
        web_env = None        
        genome_to_nuccore_response_xml = BeautifulSoup(genome_to_nuccore_response.text)
        try:
            query_key = genome_to_nuccore_response_xml.find('querykey').get_text()
            web_env = genome_to_nuccore_response_xml.find('webenv').get_text()
        except:
            if(self.verbose): print(genome_to_nuccore_response_xml)
            raise EntrezEUtilsDriverException("Query Key and/or Web Env Missing")
        if self.verbose: print("Fetching Genome fasta file from Entrez")
        get_nucleotide_sequences_response = requests.get(get_nucleotides_from_link \
            .format(**dict(query_key = query_key, web_env = web_env, form = form, mode = mode)))
        get_nucleotide_sequences_response.raise_for_status()
        # The nucleotide sequence should be located in get_genome_response.text
        if self.verbose: print("Response recieved...")
        return get_nucleotide_sequences_response.text

    def find_nucleotides_by_gene_id(self, gene_id, form = 'fasta', mode = 'text'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        if self.verbose: print("Fetching data from Entrez")
        get_nucleotide_sequences_response = requests.get(get_nucleotides_by_gene_id.format(**dict(gene_id=gene_id, form=form, mode=mode)))
        get_nucleotide_sequences_response.raise_for_status()
        # The genome sequence should be located in get_genome_response.text
        
        return get_nucleotide_sequences_response.text

    def find_gene_by_gene_id(self, gene_id, form = 'xml'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        if self.verbose: print("Fetching data from Entrez")
        get_gene_response = requests.get(get_gene_by_gene_id.format(**dict(gene_id=gene_id, form=form)))
        get_gene_response.raise_for_status()
        if self.verbose:
            open(gene_id+".gene." + form, 'w').write(get_gene_response.text)
        return get_gene_response.text

    def find_nucleotide_annotations_by_gene_id(self, gene_id, mode = 'xml'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        result = self.find_nucleotides_by_gene_id(gene_id, form = 'genbank', mode = mode)
        if self.verbose:
            open(gene_id + '.genbank.' + mode, 'w').write(result)
        return GenBankFeatureFile(result, self.opts)

    def find_protein_annotations_by_gene_id(self, gene_id, mode = 'xml'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        get_protein_sequences_response = requests.get(get_protein_by_gene_id.format(**dict(gene_id=gene_id, form='genbank', mode = mode)))
        get_protein_sequences_response.raise_for_status()
        if self.verbose:
            open(gene_id + '.peptide_record.' + mode, 'w').write(get_protein_sequences_response.text)
        return GenBankFeatureFile( get_protein_sequences_response.text, self.opts)



## GenBankFeatureFile
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the XML definition of a GenBank flat file.
class GenBankFeatureFile(object):
    def __init__(self, xml, opts):
        self.parser = BeautifulSoup(xml)
        self.opts = opts
        self.verbose = opts['verbose']
        self.chromosome = self.parser.find('iupacna').get_text()
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
        rep = "(" + ', '.join(map(repr, self.features)) + ")"
        return rep

## GenBankFeature
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the substructure of a GenBank flat file related to a single sequence feature
class GenBankFeature(object):
    def __init__(self, element, title = None):
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

        self.title = title

    def __repr__(self):
        rep = "<gi|%(gid)s (%(start)r, %(end)r %(title)s>" % self.__dict__
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
        self.end    = None
        
        self.amino_acid_sequence = element.find('iupacaa')
        # Prune later once starts and ends are set
        self.nucleotide_sequence = owner.chromosome

        self.title = element.find('seqdesc_title').get_text()
        self.annotations = map(GenBankAnnotation, element.find_all("bioseq_annot"))
        
        self._id = element.find_all('bioseq_id')
        self._desc = element.find_all("bioseq_descr")
        self._inst = element.find_all("bioseq_inst")
        self._annot = element.find_all("bioseq_annot")

    def update_genome_position(self, feature):
        self.starts = feature.starts
        self.ends   = feature.ends  
        self.start  = feature.start 
        self.end    = feature.end   

        self.nucleotide_sequence = self.nucleotide_sequence[self.start:self.end]

    def __repr__(self):
        rep = "<gi|%(gid)s Title=%(title)s, ACC=%(accession)s, ANNO=%(annotations)s>" % self.__dict__
        return rep

GenBankCDD = namedtuple("GenBankCDD", ["start", "end", "name", "definition", "e_value"])
class GenBankAnnotation(object):
    def __init__(self, element):
        self.element = element

        self.start = map(lambda x: int(x.get_text()), element.find_all("seq-interval_from"))
        self.end = map(lambda x: int(x.get_text()), element.find_all("seq-interval_to"))
        
        self.regions = [] #list(map(lambda x: x.get_text(), element.find_all("seqfeatdata_region")))
        self.name = ', '.join(map(lambda x: x.get_text(), element.find_all("prot-ref_name_e")))
        self.raw_extensions = element.find_all("seq-feat_ext")
        for raw_ext in self.raw_extensions:
            obj_type = str(raw_ext.find("object-id_str").get_text())
            obj_data = map(lambda x: x.get_text(), raw_ext.find_all("user-object_data")[0].find_all("user-field_data"))
            if obj_type == "cddScoreData":
                start = int(obj_data[0])
                end = int(obj_data[1])
                definition = str(obj_data[2]).replace('\n','')
                name = str(obj_data[3]).replace('\n','')
                e_value = float(obj_data[5])
                ext = GenBankCDD(start, end, name, definition, e_value)
                self.regions.append(ext)
            else:
                raise UnknownAnnotationException("Unknown Extension: %s" % obj_type)

    def __repr__(self):
        rep = "Annotation(Starts=%(start)s, Ends=%(end)s, Name=%(name)s, Regions=%(regions)s)" % self.__dict__
        return rep


## EntrezEUtilsDriverException
# Parent class for capturing all EUtils generated exceptions. Exception class 
# representing programmatic errors while fetching information with EntrezEUtilsDriver
class EntrezEUtilsDriverException(Exception):
    pass

## OrganismNotFoundException
# Exception indicating the organism indicated by org_name was not 
# found on Entrez
class OrganismNotFoundException(EntrezEUtilsDriverException):
    pass

class UnknownAnnotationException(EntrezEUtilsDriverException):
    pass