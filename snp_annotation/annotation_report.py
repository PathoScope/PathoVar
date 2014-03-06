import os
import re
import json
from copy import deepcopy
from collections import defaultdict

from pathovar.utils import vcf_utils, defline_parser
from pathovar.utils.fasta_utils import SequenceRecord, MutatedSequenceRecord, MutationException

def variant_to_dict(var):
    return {
        "start": var.start, "end": var.end, "ref": str(var.REF), "alts":map(str, var.ALT), 
        "quality": var.QUAL, "depth": var.INFO["DP4"],
            "_var_type": var.var_type, 
            "var_type": ("snp" if len(str(var.REF)) == len(map(str, var.ALT)[0]) 
                    else ("insertion" if len(str(var.REF)) < len(map(str, var.ALT)[0]) 
                        else "deletion"))
                }

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
            if key in self.data[org]['entries']:
                org_name = org
                break
        if org_name == None:
            raise KeyError("Key %s Not Found" % key)
        return self.data[org_name]['entries'][key]

    def __setitem__(self, key, value):
        org_name = None
        for org in self.data:
            if key in self.data[org]['entries']:
                org_name = org
                break
        if org_name == None:
            raise KeyError("Key %s Not Found" % key)
        self.data[org_name]['entries'][key] = value

    # Simplifies writing out final annotation. This forms a list of all organisms
    # being annotated. 
    def to_json_file(self):
        json_data = [deepcopy(org) for org in self.data.values()]
        for val in json_data:
            val.pop("chromosome", None)
        json.dump(json_data, open(self.vcf_path[:-3]+'json', 'w'))

    def get_annotations_from_entrez_mapping(self):
        for gene in self.genes:
            for org,val in self.annotation_dict.items():
                if gene in val.entries:
                    if val.org_name not in self.data:
                        self.data[val.org_name] = {'name':val.org_name, 
                            'accession': val.accession, 'gid': val.gid,
                            'chromosome': val.chromosome, 'entries':{}
                         }
                    if gene not in self.data:
                        entry = val.entries[gene]
                        self.data[val.org_name]['entries'][entry.gid] = entry.to_json_safe_dict()
                        self.data[val.org_name]['entries'][entry.gid]['variants'] = self.variant_by_gene[gene]

    def get_entrez_gene_annotations(self):
        if self.verbose: print("Getting Genbank Gene Comment Annotations")
        for gene in self.genes:
            gene_comments = self.annotation_manager.get_gene_comments(gene)
            self[gene]['gene_comments'] = gene_comments.to_json_safe_dict()

    def get_entrez_biosystem_pathways(self):
        if self.verbose: print("Getting Genbank BioSystem Pathway Annotations")
        for gene in self.genes:
            biosystems = self.annotation_manager.get_biosystems(gene)
            self[gene]['biosystems'] = [biosystem.to_json_safe_dict() for biosystem in biosystems]

    def generate_reference_protein_fasta_for_variants(self):
        if self.verbose: print("Generating Reference Protein Sequences.")
        sequences = []
        for org_name in self.data:
            for gene in self.data[org_name]["entries"]:
                entry = self.data[org_name]["entries"][gene]
                if entry['is_partial'] or entry['is_pseudo']:
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
                    cluster['snp_locs'] = map(variant_to_dict, current_cluster)
                    cluster['nucleotide_sequence'] = self.data[org_name]["chromosome"][(last_pos-500):last_pos]
                    cluster_mapping["intergenic-%d_%d" % (last_pos-500, last_pos)] = cluster
                    current_cluster = []
                    last_pos = var.start
                current_cluster.append(var)
            self.data[org_name]['intergenics'] = cluster_mapping

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
        for org_id, org_val in result_obj.items():
            for gene_id, gene_val in org_val.items():
                pass

    def write_text_report(self):
        handle = open(self.vcf_path[:-3] + 'report.tsv', 'w')
        for org_name, org_data in self.data.items():
            for entry_gid, entry_data in org_data['entries'].items():
                report_line = ""






