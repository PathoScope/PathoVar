import os
import re
import json

from pathovar.web.entrez_eutils import EntrezEUtilsDriver
from pathovar.web.ncbi_xml.genome_annotations import GenBankFeatureFile, CACHE_SCHEMA_VERSION as GENOME_CACHE_SCHEMA_VERSION
from pathovar.web.ncbi_xml.gene import GenbankGeneFile, CACHE_SCHEMA_VERSION as GENE_CACHE_SCHEMA_VERSION

DEFAULT_CACHE_DIR = './.anno_cache'

class EntrezAnnotationManager(object):
    def __init__(self, **opts):
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.entrez_handle = EntrezEUtilsDriver(**opts)

        # Store loaded annotations in memory
        self.genome_annotations = dict()
        self.gene_annotations = dict()

        # Try to create a data persistence cache on disk using the chosen driver
        self.cache_type = opts.get('cache_type', 'json')
        if self.cache_type == 'json':
            self.cache_manager = JSONAnnotationCacheManager(**opts)
        elif self.cache_type == None:
            self.cache_manager = None

    def get_genome_annotation(self, gid):
        data = None
        if gid in self.genome_annotations:
            data = self.genome_annotations[gid]
        elif self.cache_manager:
            try:
                data = GenBankFeatureFile( self.cache_manager.load_from_cache(gid, "genome"), verbose = self.verbose, json = True, mol_type='nucl')
                self.genome_annotations[gid] = data
            except CachedObjectMissingOrOutdatedException, e:
                data = self.find_annotations_by_gene_id(gid, xml=True, verbose = self.verbose)
                self.genome_annotations[gid] = data
                self.cache_manager.write_to_cache(data, gid, "genome")
        else: 
            data = self.find_annotations_by_gene_id(gid, xml=True, verbose = self.verbose)
            self.genome_annotations[gid] = data
        return data

    def get_gene_comments(self, gid):
        data = None
        if gid in self.gene_annotations:
            data = self.gene_annotations[gid]
        elif self.cache_manager:
            try:
                data = GenbankGeneFile(self.cache_manager.load_from_cache(gid, "gene"), verbose = self.verbose, json=True)
                self.gene_annotations[gid] = data
            except CachedObjectMissingOrOutdatedException, e:
                data = self.find_gene_comments_by_gene_id(gid, xml=True, verbose = self.verbose)
                self.gene_annotations[gid] = data
                self.cache_manager.write_to_cache(data, gid, "gene")
        else: 
            data = self.find_annotations_by_gene_id(gid, xml=True, verbose = self.verbose)
            self.gene_annotations[gid] = data
        return data

    ##
    # Merge Nucleotide and Protein file annotations together for presentation to the 
    # user. Sequence Entries between Nucleotide and Protein files are keyed in the same
    # way, but the Protein file will contain different annotations.
    def find_annotations_by_gene_id(self, gene_id, **opts):
        opts['mol_type'] = 'nucl'
        nucleotide_annotations = self.find_nucleotide_annotations_by_gene_id(gene_id, **opts)
        opts['mol_type'] = 'prot'
        # The protein genpep file for the genome does not contain region annotations
        if len(nucleotide_annotations.features) > 0:
            protein_annotations = self.find_protein_annotations_by_gene_id(nucleotide_annotations.features[0].gid, **opts)
            for gid, prot_entry in protein_annotations.entries.items():
                for annot_name, annot in prot_entry.annotations.items():
                    if annot_name in nucleotide_annotations.entries[gid].annotations:
                        nucleotide_annotations.entries[gid].annotations[annot_name].add_regions(annot.regions)
                        if nucleotide_annotations.entries[gid].annotations[annot_name].comments != (annot.comments):
                            nucleotide_annotations.entries[gid].annotations[annot_name].comments += (annot.comments)
                    else: 
                        nucleotide_annotations.entries[gid].annotations[annot_name] = annot
        else: 
            if self.verbose: print("Failed to find features for %s, data may be missing" % gene_id)
        return nucleotide_annotations

    ##
    #
    def find_nucleotide_annotations_by_gene_id(self, gene_id, **opts):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        result = self.entrez_handle.find_nucleotides_by_gene_id(gene_id, form = 'genbank', mode = 'xml')
        return GenBankFeatureFile(result, **opts)
    ##
    # 
    def find_protein_annotations_by_gene_id(self, gene_id, **opts):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        result = self.entrez_handle.find_protein_by_gene_id(gene_id, form = 'genbank', mode = 'xml')
        return GenBankFeatureFile( result, **opts)

    def find_gene_comments_by_gene_id(self, gene_id, **opts):
        gene_id = str(gene_id)
        result = self.entrez_handle.find_gene_by_gene_id(gene_id, form = 'xml')
        return GenbankGeneFile(result, **opts)

class JSONAnnotationCacheManager(object):
    GENOME_DATA = '.anno.json'
    GENE_DATA = '.gene.json'
    def __init__(self, **opts):
        self.cache_dir = opts.get("cache_dir", DEFAULT_CACHE_DIR)
        self.verbose = opts.get('verbose', False)
        self.opts = opts
        try:
            # If the cache directory is in the working directory, omit the leading os.sep
            # so that the system doesn't think it wants access to the root directory
            if self.cache_dir == os.sep + opts.get('cache_dir', '.anno_cache'):
                self.cache_dir = self.cache_dir[1:]
            os.makedirs(self.cache_dir)
        except OSError, e:
            pass

    def load_from_cache(self, query_id, data_type):
        cache_file = self.cache_dir + os.sep + 'gene_id-' + query_id
        schema_version = None
        if data_type == "genome":
            cache_file += JSONAnnotationCacheManager.GENOME_DATA
            schema_version = GENOME_CACHE_SCHEMA_VERSION
        elif data_type == "gene":
            cache_file += JSONAnnotationCacheManager.GENE_DATA
            schema_version = GENE_CACHE_SCHEMA_VERSION
        else:
            raise Exception("Annotation Data Type Missing")
        if os.path.exists(cache_file):
            json_dict = json.load(open(cache_file))
            if "schema_version" not in json_dict or (json_dict['schema_version'] != schema_version):
                if self.verbose: print("Cache Obsolete")
                raise CachedObjectMissingOrOutdatedException()
            else:
                return json_dict
        else: 
            if self.verbose: print("Cache Miss")
            raise CachedObjectMissingOrOutdatedException()

    def write_to_cache(self, data, query_id, data_type):
        cache_file = self.cache_dir + os.sep + 'gene_id-' + query_id
        if data_type == "genome":
            cache_file += JSONAnnotationCacheManager.GENOME_DATA
        elif data_type == "gene":
            cache_file += JSONAnnotationCacheManager.GENE_DATA
        else:
            raise Exception("Annotation Data Type Missing")
        json.dump(data.to_json_safe_dict(), open(cache_file, 'w'))


class CachedObjectMissingOrOutdatedException(Exception):
    pass