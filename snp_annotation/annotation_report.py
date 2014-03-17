import os
import re
import json
from copy import deepcopy
from collections import defaultdict

from pathovar.utils import vcf_utils, defline_parser
from pathovar.utils.fasta_utils import SequenceRecord, MutatedSequenceRecord, MutationException
from pathovar.web.entrez_eutils import EntrezEUtilsDriverException

GENE_DB = 'gene_comments'
BIOSYS_DB = "biosystems"

def variant_to_dict(var):
    return {
        "start": var.start, "end": var.end, "ref": str(var.REF), "alts":map(str, var.ALT), 
        "call_quality": var.QUAL, "depth": var.INFO["DP4"], "_info": var.INFO, "_format": var.FORMAT,
            "_var_type": var.var_type, 
            "var_type": ("snp" if len(str(var.REF)) == len(map(str, var.ALT)[0]) 
                    else ("insertion" if len(str(var.REF)) < len(map(str, var.ALT)[0]) 
                        else "deletion"))
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

# Perform Annotation post-processing and report building
class AnnotationReport(object):
    def __init__(self, vcf_path, annotation_manager, **opts):
        self.vcf_path = vcf_path
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.annotation_manager = annotation_manager
        self.data = dict()
        self.annotation_dict = annotation_manager.genome_annotations
        self.genes, self.variant_by_gene, self.intergenic_variants = vcf_utils.get_variant_genes(self.vcf_path)
        self.variant_by_gene = {gene:[variant_to_dict(var) for var in self.variant_by_gene[gene]] for gene in self.variant_by_gene }
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

    # Simplifies writing out final annotation. This forms a list of all organisms
    # being annotated. 
    def to_json_file(self):
        json_data = [deepcopy(org) for org in self.data.values() if len(org['entries'].keys()) > 0]
        for val in json_data:
            val.pop("chromosome", None)
        json.dump(json_data, open(self.vcf_path[:-3]+'json', 'w'))
        return self.vcf_path[:-3]+'json'

    def get_annotations_from_entrez_mapping(self):
        for org,val in self.annotation_dict.items():
            if val.org_name + ':' + val.gid not in self.data:
                self.data[val.org_name + ':' + val.gid] = {'name':val.org_name, 
                    'accession': val.accession, 'gid': val.gid,
                    'chromosome': val.chromosome, 'entries':{}
                 }
        for gene in self.genes:
            for org,val in self.annotation_dict.items():
                if gene in val.entries:
                    if gene not in self.data:
                        entry = val.entries[gene]
                        self.data[val.org_name + ':' + val.gid]['entries'][entry.gid] = entry.to_json_safe_dict()
                        self.data[val.org_name + ':' + val.gid]['entries'][entry.gid]['variants'] = self.variant_by_gene[gene]

    def get_entrez_gene_annotations(self):
        if self.verbose: print("Getting Genbank Gene Comment Annotations")
        for gene in self.genes:
            entry = self[gene]
            locus_tag = entry['gene_ref_tag']
            try:
                gene_comments = self.annotation_manager.get_gene_comments(locus_tag)
                self[gene][GENE_DB] = gene_comments.to_json_safe_dict()
            except EntrezEUtilsDriverException, e:
                if self.verbose: print(str(gene) + " is not in Gene database.")



    def get_entrez_biosystem_pathways(self):
        if self.verbose: print("Getting Genbank BioSystem Pathway Annotations")
        for gene in self.genes:
            try:
                if GENE_DB in self[gene]:
                    gene_db_id = self[gene][GENE_DB]['gene_db_id']
                    biosystems = self.annotation_manager.get_biosystems(gene_db_id)
                    self[gene][BIOSYS_DB] = [biosystem.to_json_safe_dict() for biosystem in biosystems]
                else:
                    if self.verbose: print(str(gene) + " was not queried. Missing Gene data.")
            except EntrezEUtilsDriverException, e:
                if self.verbose: print(str(gene) + " is not in BioSystem database.")


    def generate_reference_protein_fasta_for_variants(self):
        if self.verbose: print("Generating Reference Protein Sequences.")
        sequences = []
        for org_name in self.data:
            for gene in self.data[org_name]["entries"]:
                entry = self.data[org_name]["entries"][gene]
                if entry['is_partial'] or entry['is_pseudo'] or entry['is_rna']:
                    continue
                if len(entry['amino_acid_sequence']) == 0:
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

    def merge_intergenic_record_chunks(self):
        for org_name, intergenics in self.intergenic_variants.items():
            cluster_mapping = dict()
            current_cluster = []
            last_pos = intergenics[0].start
            for var in intergenics:
                if var.start - 500 >= last_pos:
                    cluster = dict()
                    cluster['start'] = current_cluster[0].start - 250
                    cluster['upstream_id'] = current_cluster[0].INFO['GENE'].split(":")[1].split("~")[0]
                    cluster['downstream_id'] = current_cluster[-1].INFO['GENE'].split(":")[1].split("~")[0]
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
                    current_cluster = []
                    last_pos = var.start
                current_cluster.append(var)
            self.data[org_name]['entries'].update(cluster_mapping)


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



    ## TODO
    ## Mutation transformation validation fails on sequences with indels, and indices are unreliable
    ## so it is not included in the pipeline
    def generate_mutant_nucleotide_sequences(self):
        if self.verbose: print("Generating Mutant Nucleotide Sequences.")
        sequences = []
        errs = defaultdict(list)
        for org_name in self.data:
            if self.verbose: print("Handling %s" % org_name)
            for gene in self.data[org_name]["entries"].values():
                    defline = 'gi|%(gid)s|ref|%(accession)s| %(title)s' % gene
                    #if self.verbose: print("\tHandling %s" % defline)
                    seq_rec = MutatedSequenceRecord(defline, gene['nucleotide_sequence'], defline_parser, **gene)
                    if gene['is_partial'] or gene['is_pseudo']:
                        continue
                    if len(gene['nucleotide_sequence']) == 0:
                        if self.verbose: print("Skipping %s due to missing nucleotide data" % gene['gid'])
                        seq_rec.defline += "_" + org_name
                        errs[org_name].append(seq_rec)
                        continue
                    try:
                        seq_rec.sweep_mutations()
                        gene['mutant_nucleotide_sequence'] = seq_rec.mutated_sequence
                        gene['mutant_nucleotide_sequence_variant_indices'] = seq_rec.mutated_indices
                        sequences.append(seq_rec)
                    except MutationException, e:
                        seq_rec.defline += "_" + org_name
                        errs[org_name].append((seq_rec, e,))
                        if self.verbose: print(seq_rec, e)
                            
        fasta_name = self.vcf_path[:-3] + "variant_mutant_nucleotides.fa"
        fasta_handle = open(fasta_name, 'w')
        for seq in sequences:
            fasta_handle.write(seq.to_fasta_format())
        fasta_handle.close()

        if len(errs.values()) > 0:
            with open(self.vcf_path[:-3] + "err", 'w') as err_log:
                for org_err in errs.values():
                    for seq_err in org_err:
                        err_log.write(','.join(map(str, seq_err)) + '\n')
        return fasta_name

    def consume_blast_results(self,db_name,blast_results):
        if self.verbose: print("Consuming %s BLAST Results." % db_name)
        for query in blast_results.queries:
            query_entry = self[query]
            if not 'blast_hits' in query_entry:
                query_entry['blast_hits'] = dict()
            query_entry['blast_hits'][db_name] = blast_results.queries[query].to_json_safe_dict()['hits']

    def consume_generic_result(self, result_obj):
        if self.verbose: print("Consuming %s Results" % result_obj.name)
        for org_name, org_val in result_obj.items():
            for gene_name, gene_val in org_val.items():
                pass






