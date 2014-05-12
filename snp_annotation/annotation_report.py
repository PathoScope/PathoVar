import os
import re
import json
from copy import deepcopy
from collections import defaultdict

from pathovar.utils import vcf_utils, defline_parser
from pathovar.utils.fasta_utils import SequenceRecord, MutatedSequenceRecord, MutationException
from pathovar.web.entrez_eutils import EntrezEUtilsDriverException
from pathovar.snp_caller import compute_sam_coverage

GENE_DB = 'gene_comments'
BIOSYS_DB = "biosystems"

def variant_to_dict(var):
    return {
        "start": var.start, "end": var.end, "ref": str(var.REF), "alts":map(str, var.ALT), 
        "call_quality": var.QUAL, "depth": var.INFO["DP4"], "_info": var.INFO, "_format": var.FORMAT,
            "_var_type": var.var_type, 
            "var_type": ("snp" if len(str(var.REF)) == len(map(str, var.ALT)[0]) 
                    else ("insertion" if len(str(var.REF)) < len(map(str, var.ALT)[0]) 
                        else "deletion")),
            "eff": {
                "type": "UNKNOWN" if len(str(var.REF)) == len(map(str, var.ALT)[0]) else "FRAME_SHIFT"
            }
        }

# Translate an entry to it's most common identifiers
def gene_name(entry):
    if entry['gene_symbol'] is not None:
        return entry['gene_symbol']
    elif entry['gene_ref_tag'] is not None:
        return entry['gene_ref_tag']
    else:
        return entry['gid']

def spans_variant(entry, start, stop):
    vars_spanned = []
    for i, variant in enumerate(entry['variants']):
        for start_pos in start:
            for stop_pos in stop:
                if ((variant['start'] - entry['start'] <= start_pos and variant['end'] - entry['start'] >= stop_pos) or
                (variant['start'] - entry['start'] >= start_pos and variant['end'] - entry['start'] <= stop_pos) or
                (variant['start'] - entry['start'] <= start_pos and variant['end'] - entry['start'] >= start_pos) or
                (variant['start'] - entry['start'] <= stop_pos and variant['end'] - entry['start'] >= stop_pos)):
                    vars_spanned.append(i)

    return False if len(vars_spanned) == 0 else vars_spanned

def score_heuristic(entry, snp_score = .1, missense_score = 2, frame_shift_score = 2/3.0, min_quality = 20, 
                    blast_hit_score = 2, multi_drug_bonus = 2, uncov_region_score = 2):
    variant_score = 1
    seq_len = float(entry['end'] - entry['start'])
    for variant in entry["variants"]:
        if variant['call_quality'] < min_quality:
            continue
        if variant["eff"]["type"] in ["UNKNOWN", "SYNONYMOUS_CODING", "INTERGENIC"]:
            variant_score += snp_score

        elif variant["eff"]["type"] in ["NON_SYNONYMOUS_CODING", "NON_SYNONYMOUS_START", \
                                        "SYNONYMOUS_STOP", "NON_SYNONYMOUS_STOP"]:
            variant_score += missense_score

        elif variant["eff"]["type"] in ["STOP_GAINED", "STOP_LOST", "START_LOST", "RARE_AMINO_ACID"]:
            variant_score += seq_len/3.0 * min(seq_len / 2000, 1)

        elif variant["eff"]["type"] in ["FRAME_SHIFT", "SPLICE_SITE_ACCEPTOR", "SPLICE_SITE_DONOR"]:
            start_point = variant['start'] - entry['start']
            percent_shift = (start_point / seq_len)
            # The earlier the shift occurs, the greater the score multiplier, 
            # but penalize short sequences ( < 1000 bp) 
            score = (frame_shift_score / percent_shift) * min(seq_len / 2000, 1)
            variant_score += score
            #variant_score += frame_shift_score
        else:
            print("Variant Effect Not Recognized: %s" % variant["eff"]["type"])
    
    blast_score = 1
    for blast_db, blast_hits in entry['blast_hits'].items():
        for hit in blast_hits:
            blast_score += blast_hit_score
            if re.search(r'multidrug', hit['hit_def'], re.IGNORECASE):
                blast_score += multi_drug_bonus

    coverage_score = 1
    for region in entry["uncovered_regions"]:
        # if(entry['mean_coverage'] == 0):
        #     uncov_region_score * --
        coverage_score += uncov_region_score

    entry["score"] = variant_score * blast_score * coverage_score
    return entry

def gene_filter(entry, hypothetical = True, intergenic = True, min_quality = 20):
    if re.search(r"hypothetical", entry['title']) and hypothetical:
        return False
    elif entry['is_intergenic'] and intergenic:
        return False
    if any([var['call_quality'] < min_quality for var in entry['variants'] if min_quality is not None]):
        return False

    return True



# Makes sure that all entries have all of the possible "extra"
# record data included. 
def normalize_entry_model(entry):
    if "mean_coverage" not in entry:
        entry["mean_coverage"] = None
    if "uncovered_regions" not in entry:
        entry["uncovered_regions"] = []
    if "is_intergenic" not in entry:
        entry["is_intergenic"] = False
    if "is_partial" not in entry:
        entry["is_partial"] = False
    if "is_pseudo" not in entry:
        entry["is_pseudo"] = False
    if "variants" not in entry:
        entry["variants"] = []
    if "amino_acid_sequence" not in entry:
        entry["amino_acid_sequence"] = ""
    if "blast_hits" not in entry:
        entry["blast_hits"] = {}
    for variant in entry['variants']:
        if "eff" not in variant:
            variant['eff']={"type":"UNKNOWN"}


# Perform Annotation post-processing and report building. Portions of this class
# are written assuming it only had access to the annotated VCF file generated by 
# the VariantLocator class. Other sections assume it has access to the data objects
# used to create the annotated VCF. There is no reason not to assume you have access
# to the data objects. 
class AnnotationReport(object):
    def __init__(self, vcf_path = None, variant_locator = None, annotation_manager = None, **opts):
        if(annotation_manager is None):
            raise AnnotationReportException("AnnotationReport instance requires an AnnotationManager instance, did not recieve one.")
        if(vcf_path is None and variant_locator is None):
            raise AnnotationReportException("AnnotationReport instance requires an annotated VCF file path and a VariantLocator instance, but recieved neither.")
        self.variant_locator = variant_locator
        if(self.variant_locator is not None):
            self.vcf_path = variant_locator.vcf_file
        else:
            self.vcf_path = vcf_path
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.annotation_manager = annotation_manager
        self.data = dict()
        self.annotation_dict = annotation_manager.genome_annotations
        # Parse annotated VCF from previous scan. This stage could be instead performed on the data objects from
        # an instance of VariantLocator
        if(self.variant_locator is not None):
            self.genes, self.variant_by_gene, self.intergenic_variants = variant_locator.get_variant_genes()
        elif(self.vcf_path is not None):
            self.genes, self.variant_by_gene, self.intergenic_variants = vcf_utils.get_variant_genes(self.vcf_path)
        else:
            raise AnnotationReportException("No annotated variant information provided.")
        # Transform the variants per gene dictionary from vcf._Record to dictionaries to serialize as JSON
        self.variant_by_gene = {gene:[variant_to_dict(var) for var in self.variant_by_gene[gene]] for gene in self.variant_by_gene }
        # Copy data from annotation_manager to dictionary format to freely decorate with new annotations and serialize as JSON
        self.get_annotations_from_entrez_mapping()

    # Abstraction that hides the fact we are dealing with multiple organisms
    def __getitem__(self, key):
        org_name = None
        for org in self.data:
            #print(org, (self.data[org]['entries']))
            if key in self.data[org]['entries']:
                org_name = org
                break
        if org_name == None:
            raise KeyError(key)
        return self.data[org_name]['entries'][key]

    def __setitem__(self, key, value):
        org_name = None
        for org in self.data:
            if key in self.data[org]['entries']:
                org_name = org
                break
        if org_name == None:
            raise KeyError(key)
        self.data[org_name]['entries'][key] = value

    def __iter__(self):
        for gid in self.genes:
            yield self[gid]

    def get_org_by_gid(self, gid):
        for org_name, org in self.data.items():
            if org["gid"] == gid:
                return org_name
        raise KeyError("No Organism/Genome with GID:%s" % gid)

    # Simplifies writing out final annotation. This forms a list of all organisms
    # being annotated. 
    def to_json_file(self):
        if self.verbose: print("Saving results")
        json_data = [deepcopy(org) for org in self.data.values() if len(org['entries'].keys()) > 0]
        for organism in json_data:
            for gid, entry in organism['entries'].items():
                normalize_entry_model(entry)
                
        for val in json_data:
            val.pop("chromosome", None)
        json.dump(json_data, open(self.vcf_path[:-3]+'json', 'w'))
        return self.vcf_path[:-3]+'json'

    ##
    # Copy all of the relevant data from the Entrez data in the AnnotationManager
    # into JSON safe dictionary format for attaching new annotations to. 
    def get_annotations_from_entrez_mapping(self):
        for org,val in self.annotation_dict.items():
            if val.org_name + ':' + val.gid not in self.data:
                self.data[val.org_name + ':' + val.gid] = {
                    'name':val.org_name, "sub_name": val.sub_name,
                    'accession': val.accession, 'gid': val.gid,
                    'db_tag_data': val.db_tag_data,
                    'chromosome': val.chromosome, 'entries':{}
                 }
        for gene in self.genes:
            for org,val in self.annotation_dict.items():
                if gene in val.entries:
                    if gene not in self.data:
                        entry = val.entries[gene]
                        self.data[val.org_name + ':' + val.gid]['entries'][entry.gid] = entry.to_json_safe_dict()
                        self.data[val.org_name + ':' + val.gid]['entries'][entry.gid]['variants'] = self.variant_by_gene[gene]

    ##
    # Get Entrez Gene annotations by locus tag
    def get_entrez_gene_annotations(self):
        if self.verbose: print("Getting Genbank Gene Comment Annotations")
        for gene in self.genes:
            entry = self[gene]
            if entry['is_intergenic']:
                continue
            if 'gene_ref_tag' not in entry:
                entry['gene_ref_tag'] = None
            locus_tag = entry['gene_ref_tag']
            if locus_tag is None:
                #if self.verbose: print(str(gene) + " is not in Gene database (No Locus Tag)")
                continue
            try:
                gene_comments = self.annotation_manager.get_gene_comments(locus_tag)
                self[gene][GENE_DB] = gene_comments.to_json_safe_dict()
            except EntrezEUtilsDriverException, e:
                pass
                #if self.verbose: print(str(gene) + " is not in Gene database. (No Record)")

    ##
    # Get Entrez BioSystems (like KEGG Pathways) using the Entrez Gene ID of the sequence.
    # If the sequence does not have such an ID, then skip it, it hasn't been well studied enough
    # to have a BioSystem entry.
    def get_entrez_biosystem_pathways(self):
        if self.verbose: print("Getting Genbank BioSystem Pathway Annotations")
        for gene in self.genes:
            try:
                if GENE_DB in self[gene]:
                    gene_db_id = self[gene][GENE_DB]['gene_db_id']
                    biosystems = self.annotation_manager.get_biosystems(gene_db_id)
                    self[gene][BIOSYS_DB] = [biosystem.to_json_safe_dict() for biosystem in biosystems]
                #else:
                #     if self.verbose: print(str(gene) + " was not queried. Missing Gene data.")
            except EntrezEUtilsDriverException, e:
                if self.verbose: print(str(gene) + " is not in BioSystem database. (No Record)")


    ##
    # For each gene that contains one or more variants, 
    def generate_reference_protein_fasta_for_variants(self):
        if self.verbose: print("Generating Reference Protein Sequences.")
        sequences = []
        for org_name in self.data:
            for gene in self.data[org_name]["entries"]:
                entry = self.data[org_name]["entries"][gene]
                if entry['is_partial'] or entry['is_pseudo'] or entry['is_rna'] or entry['is_intergenic']:
                    continue
                if entry['amino_acid_sequence'] is None or len(entry['amino_acid_sequence']) == 0:
                    if self.verbose: print(str(entry['gid']) + " Protein Sequence Missing")
                    continue
                defline = 'gi|%(gid)s|ref|%(accession)s| %(title)s' % entry
                seq_rec = SequenceRecord(defline, entry['amino_acid_sequence'][0], defline_parser)
                sequences.append(seq_rec)
        fasta_name = self.vcf_path[:-3] + "variant_ref_protein.fa"
        fasta_handle = open(fasta_name, 'w')
        for seq in sequences:
            fasta_handle.write(seq.to_fasta_format())
        fasta_handle.close()
        return fasta_name


    ##
    # Collect all intergenic variants and assemble them into 500 nt spans
    # between genes. 
    def merge_intergenic_record_chunks(self):
        for org_tag, intergenics in self.intergenic_variants.items():
            # Interpreting the VCF File means loss of data. Use the gid after the colon
            # in the encoded string to get a key to get the full organism name
            org_gid = org_tag.split(":")[-1]
            org_name = self.get_org_by_gid(org_gid)
            cluster_mapping = dict()
            current_cluster = []
            last_pos = intergenics[0].start
            for var in intergenics:
                if var.start - 500 >= last_pos:
                    cluster = dict()
                    cluster['start'] = current_cluster[0].start - 250
                    # If data source is a dictionary, the data must be parsed from the info line
                    if(type(current_cluster[0].INFO["GENE"]) == dict):
                        cluster['upstream_id'] = current_cluster[0].INFO['GENE'].split(":")[1].split("~")[0]
                        cluster['downstream_id'] = current_cluster[-1].INFO['GENE'].split(":")[1].split("~")[0]
                    # Else the data source is an object that does not need to be parsed
                    else:
                        cluster['upstream_id'] = current_cluster[0].annotations[0].upstream_id
                        cluster['downstream_id'] = current_cluster[0].annotations[0].upstream_id
                    cluster['end'] = current_cluster[0].start + 750
                    cluster['variants'] = map(variant_to_dict, current_cluster)
                    cluster['nucleotide_sequence'] = self.data[org_name]["chromosome"][(cluster['start']):cluster['end']]
                    cluster['is_intergenic'] = True
                    cluster['is_rna'] = False
                    cluster['is_partial'] = False
                    cluster['is_pseudo'] = False
                    cluster['title'] = "%s-intergenic-%d_%d-%s" % (cluster['upstream_id'], (cluster['start']), cluster['start'], 
                        cluster['downstream_id'])
                    cluster['gid'] = self.data[org_name]['gid']
                    cluster['accession'] = self.data[org_name]['accession']
                    cluster_mapping[cluster['title']] = cluster
                    self.genes[cluster['title']] = {'gi': cluster['title'], 
                                                    'acc': self.data[org_name]['accession'],
                                                    'strand': '?',
                                                    'title': cluster['title']
                                                    }
                    current_cluster = []
                    last_pos = var.start
                current_cluster.append(var)
            self.data[org_name]['entries'].update(cluster_mapping)


    ##
    # For each gene, compute for each annotation whether it spans any of the variants found
    # for its gene. 
    def compute_regions_spanning_variants(self):
        for org_name in self.data:
            for gene in self.data[org_name]['entries']:
                entry = self.data[org_name]['entries'][gene]
                try:
                    for annot in entry["annotations"].values():
                        start_positions = annot["start"] if len(annot['start']) > 0 else annot["regions"][0]['from']
                        end_positions = annot["end"] if len(annot['end']) > 0 else annot["regions"][0]['to']
                        start_positions = map(lambda x: x * 3, start_positions)
                        end_positions = map(lambda x: x * 3, end_positions)
                        spanned = spans_variant(entry, start_positions, end_positions)
                        annot['spans_variant'] = spanned
                except KeyError, e:
                    # Not all Entries have an annotations field.
                    pass

    ##
    # Consume the JSON serialization from a BlastResultParser and associate
    # each query gene entry with its hits. 
    def consume_blast_results(self,db_name,blast_results):
        if self.verbose: print("Consuming %s BLAST Results." % db_name)
        for query in blast_results.queries:
            query_entry = self[query]
            if not 'blast_hits' in query_entry:
                query_entry['blast_hits'] = dict()
            query_entry['blast_hits'][db_name] = blast_results.queries[query].to_json_safe_dict()['hits']


    ##
    # Consume snpEff results and update each gene's variant dictionaries to contain
    # effect data where available.
    def consume_snpeff_results(self, snpeff_results):
        for genome, effects in snpeff_results.items():
            for locus, locus_effect in effects.items():
                if 'gID' in locus_effect:
                    gene = self[locus_effect['gID']]
                    for nth, nth_effect in locus_effect.items():
                        if nth == 'gID': continue
                        # VCF POS is 1-indexed
                        variants = [var for var in gene['variants'] if (var['start'] == (int(locus) - 1))]
                        if len(variants) == 0:
                            raise Exception("No matching variant found when mapping snpEff results")
                        var = variants[0]
                        var['eff'] = nth_effect


    def normalize_all_entries(self):
        if self.verbose : print("Normalizing all entry structures")
        for org_name, org_data in self.data.items():
            for gid, entry in org_data['entries'].items():
                normalize_entry_model(entry)

    def score_all_entries(self):
        if self.verbose : print("Running scoring heuristic")
        for org_name, org_data in self.data.items():
            for gid, entry in org_data['entries'].items():
                score_heuristic(entry)

    ##
    # Consume coverage data to get mean coverage for each gene that contains a variant
    # and place regions without coverage within genes or intergenic regions. Creates new 
    # entries when an uncovered region has not previously been captured. 
    def compute_coverage_span(self, coverage_data):
        if self.verbose:
            print("Computing base coverage")
        coverage_dict = coverage_data['coverage_dict']
        uncovered_regions_dict = coverage_data['uncovered_regions_dict']
        for org_name, org_val in self.data.items():
            org_gid = org_val['gid']
            if org_gid not in coverage_dict: 
                continue
            for gene, entry in org_val['entries'].items():
                entry['uncovered_regions'] = []
                try:
                    # Not all entries will have start/end points
                    coverage_score = compute_sam_coverage.extract_mean_coverage(coverage_dict[org_gid], entry['start'], entry['end'])
                    entry['mean_coverage'] = coverage_score
                except KeyError, e:
                    pass
            genome_annotations = self.annotation_manager.get_genome_annotation(org_gid)
            intergenic_count = 0
            current_cluster = []
            for region in uncovered_regions_dict[org_gid]:
                locations = genome_annotations.locate_snp_site(region, forward_check = True)
                #print("Forward search of uncovered regions", len(locations))
                #region.location = [loc.gid for loc in locations]
                for location in locations:
                    ##
                    # This doesn't accomodate placing uncovered regions in already known variant intergenic
                    # regions. There are also so many of them that it is not worth the dead weight
                    if location.is_intergenic:
                        intergenic_count += 1
                        # # Build new intergenic record
                        # entry = dict()
                        # entry['start'] = region.start - 50
                        # entry['end'] =  region.end + 50
                        # entry['gid'] = org_gid
                        # entry['accession'] = org_val['accession']
                        # entry['upstream_id'] = location.upstream_id
                        # entry['downstream_id'] = location.downstream_id
                        # entry['is_intergenic'] = True
                        # entry['is_rna'] = False
                        # entry['is_partial'] = False
                        # entry['is_pseudo'] = False
                        # entry['nucleotide_sequence'] = self.data[org_name]["chromosome"][(entry['start']):entry['end']]
                        # entry['title'] = "%s-intergenic-%d_%d-%s" % (entry['upstream_id'], (entry['start']), entry['start'], 
                        #                 entry['downstream_id'])
                        # entry['variants'] = []
                        # entry['uncovered_regions']=[region.to_json_safe_dict()]

                        # org_val['entries'][entry['title']] = entry
                    # In an existing gene entry
                    elif location.gid in org_val['entries']:
                        entry = self[location.gid]
                        if 'uncovered_regions' not in entry:
                            entry['uncovered_regions'] = []
                        entry['uncovered_regions'].append(region.to_json_safe_dict())
                    else:
                        # Build new gene entry
                        org_val['entries'][location.gid] = location.to_json_safe_dict()
                        self.genes[location.gid] = {
                            "acc": location.accession,
                            "gi": location.gid,
                            "strand": location.strand,
                            "title": location.title
                        }
                        entry = org_val['entries'][location.gid]
                        if 'uncovered_regions' not in org_val['entries'][location.gid]:
                            entry['uncovered_regions'] = []

                        entry['uncovered_regions'].append(region.to_json_safe_dict())
                        coverage_score = compute_sam_coverage.extract_mean_coverage(coverage_dict[org_gid], entry['start'], entry['end'])
                        entry['mean_coverage'] = coverage_score
                        entry['variants'] = []
                

    ##
    # TODO Incomplete. 
    # Consume a nested dictionary of results and associate the results for a
    # given gene with its entry. 
    def consume_generic_result(self, result_obj):
        if self.verbose: print("Consuming %s Results" % result_obj.name)
        for org_name, org_val in result_obj.items():
            for gene_name, gene_val in org_val.items():
                pass

class AnnotationReportException(Exception):
    pass






