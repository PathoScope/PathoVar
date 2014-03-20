import os
import re
import json

from pathovar.web.entrez_eutils import EntrezEUtilsDriver
from pathovar.web.ncbi_xml.genome_annotations import GenBankFeatureFile
from pathovar.web.ncbi_xml.gene import GenBankGeneFile, GenBankBioSystemFile, GenBankGeneToBioSystem

class EntrezAnnotationManager(object):
    def __init__(self, **opts):
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.entrez_handle = EntrezEUtilsDriver(**opts)

        # Store loaded annotations in memory
        self.genome_annotations = dict()
        self.gene_annotations = dict()
        self.biosystem_annotations = dict()

        # Try to create a data persistence cache on disk using the chosen driver
        self.cache_type = opts.get('cache_type', 'json')
        if self.cache_type == 'json':
            self.cache_manager = JSONAnnotationCacheManager(**opts)
        elif self.cache_type == None:
            self.cache_manager = None

    def genome_to_accesion_and_codon_table(self):
        mapping = {gid:{'accession': genome.accession, "codon_table":genome.genetic_code} for gid, genome in self.genome_annotations.items()}
        return mapping

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

    def get_gene_comments(self, locus_tag):
        data = None
        if locus_tag in self.gene_annotations:
            data = self.gene_annotations[locus_tag]
        elif self.cache_manager:
            try:
                data = GenBankGeneFile(self.cache_manager.load_from_cache(locus_tag, "gene"), verbose = self.verbose, json=True)
                self.gene_annotations[locus_tag] = data
            except CachedObjectMissingOrOutdatedException, e:
                data = self.find_gene_comments_by_locus_tag(locus_tag, xml=True, verbose = self.verbose)
                self.gene_annotations[locus_tag] = data
                self.cache_manager.write_to_cache(data, locus_tag, "gene")
        else: 
            data = self.find_annotations_by_gene_id(gid, xml=True, verbose = self.verbose)
            self.gene_annotations[gid] = data
        return data

    def get_biosystems(self, gene_db_id):
        biosystem_ids = None
        if self.cache_manager:
            try:
                gene_to_biosystem = GenBankGeneToBioSystem(self.cache_manager.load_from_cache(gene_db_id, "gene-biosystem"), json=True)
                biosystem_ids = gene_to_biosystem.biosystem_ids
            except CachedObjectMissingOrOutdatedException:
                gene_to_biosystem = GenBankGeneToBioSystem(self.entrez_handle.find_biosystem_ids_by_gene_db_id(gene_db_id), gene_id = gene_db_id, xml = True)
                self.cache_manager.write_to_cache(gene_to_biosystem, gene_db_id, "gene-biosystem")
                biosystem_ids = gene_to_biosystem.biosystem_ids
        else:
            gene_to_biosystem = GenBankGeneToBioSystem(self.entrez_handle.find_biosystem_ids_by_gene_db_id(gene_db_id), gene_id = gene_db_id, xml = True)
            biosystem_ids = gene_to_biosystem.biosystem_ids
        biosystems = []
        for bsid in biosystem_ids:
            data = None
            if bsid in self.biosystem_annotations:
                data = self.biosystem_annotations[bsid]
            elif self.cache_manager:
                try:
                    data = GenBankBioSystemFile(self.cache_manager.load_from_cache(bsid, "biosystem"), verbose = self.verbose, json=True)
                    self.biosystem_annotations[bsid] = data
                except CachedObjectMissingOrOutdatedException, e:
                    data = GenBankBioSystemFile(self.entrez_handle.find_biosystem_by_bsid(bsid), xml=True, verbose = self.verbose)
                    self.biosystem_annotations[bsid] = data
                    self.cache_manager.write_to_cache(data, bsid, "biosystem")
            else:
                data = GenBankBioSystemFile(self.entrez_handle.find_biosystem_by_bsid(bsid), xml=True, verbose = self.verbose)
                self.biosystem_annotations[bsid] = data
            biosystems.append(data)
        return biosystems

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

    def find_gene_comments_by_locus_tag(self, locus_tag, **opts):
        result = self.entrez_handle.find_gene_by_locus_tag(locus_tag)
        return GenBankGeneFile(result, **opts)

class AnnotationCacheManagerBase(object):
    DEFAULT_CACHE_DIR = "./.anno_cache"
    def __init__(self, **opts):
        self.cache_dir = opts.get("cache_dir", AnnotationCacheManagerBase.DEFAULT_CACHE_DIR)
        self.verbose = opts.get('verbose', False)


class JSONAnnotationCacheManager(object):
    GENOME_DATA = '.anno.json'
    GENE_DATA = '.gene.json'
    GENE_BIOSYSTEM_DATA = '.gene-biosys.json'
    BIOSYSTEM_DATA = '.biosys.json'

    def __init__(self, **opts):
        self.cache_dir = opts.get("cache_dir", AnnotationCacheManagerBase.DEFAULT_CACHE_DIR)
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
            schema_version = GenBankFeatureFile.CACHE_SCHEMA_VERSION
        elif data_type == "gene":
            cache_file += JSONAnnotationCacheManager.GENE_DATA
            schema_version = GenBankGeneFile.CACHE_SCHEMA_VERSION
        elif data_type == "biosystem":
            cache_file += JSONAnnotationCacheManager.BIOSYSTEM_DATA
            schema_version = GenBankBioSystemFile.CACHE_SCHEMA_VERSION
        elif data_type == "gene-biosystem":
            cache_file += JSONAnnotationCacheManager.GENE_BIOSYSTEM_DATA
            schema_version = GenBankGeneToBioSystem.CACHE_SCHEMA_VERSION
        else:
            raise Exception("Annotation Data Type Missing")
        if os.path.exists(cache_file):
            json_dict = None
            try:
                json_dict = json.load(open(cache_file))
            except ValueError, e:
                os.remove(cache_file)
                raise CachedObjectMissingOrOutdatedException()
            if "schema_version" not in json_dict or (json_dict['schema_version'] != schema_version):
                #if self.verbose: print("Cache Obsolete")
                raise CachedObjectMissingOrOutdatedException()
            else:
                return json_dict
        else: 
            #if self.verbose: print("Cache Miss")
            raise CachedObjectMissingOrOutdatedException()

    def write_to_cache(self, data, query_id, data_type):
        cache_file = self.cache_dir + os.sep + 'gene_id-' + query_id
        if data_type == "genome":
            cache_file += JSONAnnotationCacheManager.GENOME_DATA
        elif data_type == "gene":
            cache_file += JSONAnnotationCacheManager.GENE_DATA
        elif data_type == "biosystem":
            cache_file += JSONAnnotationCacheManager.BIOSYSTEM_DATA
        elif data_type == "gene-biosystem":
            cache_file += JSONAnnotationCacheManager.GENE_BIOSYSTEM_DATA
        else:
            raise Exception("Annotation Data Type Missing")
        try:
            json.dump(data.to_json_safe_dict(), open(cache_file, 'w'))
        except TypeError, e:
            try:
                if(self.verbose): print("Cache Write Error", e)
                os.remove(cache_file)
            except OSError, os_err:
                if(self.verbose): print(os_err)


class CachedObjectMissingOrOutdatedException(Exception):
    pass