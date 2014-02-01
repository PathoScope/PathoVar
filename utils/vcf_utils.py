import re
import os
import json
from collections import defaultdict
# 3rd Party Imports
import vcf
from vcf.parser import _Filter
from vcf.filters import Base as VCFFilterBase

from pathovar.utils import defline_parser
from pathovar.utils.fasta_utils import SequenceRecord

## filter_vcf
# Based on vcf_filter.py from PyVCF. Compatible with instantiated
# PyVCF _Filter objects.
# @param keep If True, keep sequences that pass a filter, else exclude sequences that pass the filter
def filter_vcf(file_path, filters, keep = True, short_circuit = False):
    inp = vcf.Reader(open(file_path, 'r'))

    # build filter chain
    chain = []
    for filter_obj in filters:
        chain.append(filter_obj)
        short_doc = filter_obj.__doc__ or ''
        short_doc = short_doc.split('\n')[0].lstrip()
        # add a filter record to the output
        inp.filters[filter_obj.filter_name()] = _Filter(filter_obj.filter_name(), short_doc)

    # output must be created after all the filter records have been added
    output = vcf.Writer(open(file_path + '.filt.vcf', 'w'), inp)

    # apply filters
    for record in inp:
        output_record = []
        for filt in chain:
            result = filt(record)
            if result is None:
                continue
            
            output_record.append(bool(result))
            record.add_filter(filt.filter_name())

        output_record = all(output_record) and len(output_record) > 0
        if not keep:
            output_record = not output_record
        if output_record:
            # use PASS only if other filter names appear in the FILTER column
            if record.FILTER is None and not keep: record.FILTER = 'PASS'
            output.write_record(record)

def filter_vcf_in_memory(variant_reader, filters, keep=True, short_circuit = False):
    # build filter chain
    chain = []
    for filter_obj in filters:
        chain.append(filter_obj)
        short_doc = filter_obj.__doc__ or ''
        short_doc = short_doc.split('\n')[0].lstrip()
        # add a filter record to the output
        variant_reader.filters[filter_obj.filter_name()] = _Filter(filter_obj.filter_name(), short_doc)

    keepers = []
    # apply filters
    for record in variant_reader:
        output_record = []
        for filt in chain:
            result = filt(record)
            
            output_record.append(bool(result))
            record.add_filter(filt.filter_name())

        output_record = all(output_record) and len(output_record) > 0
        if not keep:
            output_record = not output_record
        if output_record:
            keepers.append(record)

    return keepers

def get_variant_genes(vcf_path):
    genes = {}
    variant_by_gene = defaultdict(list)
    reader = vcf.Reader(open(vcf_path,'r'))
    for var in reader:
        gene_info = var.INFO.get('GENE', '')
        if gene_info == '': continue
        
        gene_info_groups = gene_info.split('||')
        general_info = gene_info_groups[0].split('|')
        
        # Fields are key:value pairs. Certain fields have special handling
        # The first position will have a leading paren to exclude
        general_info[0] = general_info[0][1:]
        general_info = [pair.split(':') for pair in general_info]
        gene_dict = {val[0]: val[1] for val in general_info}
        for key in gene_dict:
            if key != 'acc':
                gene_dict[key] = gene_dict[key].replace('_', ' ')
        genes[gene_dict['gi']] = gene_dict
        variant_by_gene[gene_dict['gi']].append(var)

    return(genes, variant_by_gene)

def vcf_to_gene_report(vcf_path):
    genes, variant_by_gene = get_variant_genes(vcf_path)

    report = open(os.path.splitext(vcf_path)[0] + '.variant_report.tsv', 'w')
    report.write('Gene\tVariant Positions\n')
    for gene in genes.values():
        line = "gi|{gi}|" 
        if 'acc' in gene:
            line +="ref|{acc}|" 
        line += " {title}"
        line = line.format(**gene)
        for var in variant_by_gene[gene['gi']]:
            line += '\t({POS}:{REF}->{ALT})'.format(**var.__dict__)
        report.write(line + '\n')
    report.close()

def generate_html_report_json(vcf_path, annotation_dict):
    genes, variant_by_gene = get_variant_genes(vcf_path)
    variant_by_gene = {gene:[{"start": var.start, "end": var.end, "ref": str(var.REF), "alts":map(str, var.ALT)} for var in variant_by_gene[gene]] for gene in variant_by_gene}
    json_dict = {}
    for gene in genes:
        for org,val in annotation_dict.items():
            if gene in val.entries:
                if val.org_name not in json_dict:
                    json_dict[val.org_name] = {'name':val.org_name,'entries':{}}
                if gene not in json_dict:
                    entry = val.entries[gene]
                    json_dict[val.org_name]['entries'][entry.gid] = entry.to_json_safe_dict()
                    json_dict[val.org_name]['entries'][entry.gid]['variants'] = variant_by_gene[gene]
    # Drop the keys, passing only a list of unique reference strain annotations
    json.dump(json_dict.values(), open(vcf_path[:-3] + 'json', 'w'))
    return vcf_path[:-3] + 'json'


def generate_reference_fasta_for_variants(vcf_path, annotation_dict):
    genes, variant_by_gene = get_variant_genes(vcf_path)

    sequences = []
    for gene in genes:
        for val in annotation_dict.values():
            if gene in val.entries:
                entry = val.entries[gene]
                assert entry.nucleotide_sequence != ""
                defline = 'gi|%(gid)s|ref|%(accession)s| %(title)s' % entry.__dict__
                seq_rec = SequenceRecord(defline, entry.nucleotide_sequence, defline_parser)
                sequences.append(seq_rec)
    fasta_name = vcf_path[:-3] + "variant_refs.fa"
    fasta_handle = open(fasta_name, 'w')
    for seq in sequences:
        fasta_handle.write(seq.to_fasta_format())
    fasta_handle.close()
    return fasta_name




## VCF Filter Classes
class FilterByChromMatch(VCFFilterBase):
    '''Filter a VCF File by regular expresion match over its CHROM column'''
    name = "regmatch"
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--pattern', type=str, default='.*',
                help='Regular Expression to match the CHROM column')

    def __init__(self, args):
        self.pattern = args.pattern

    def __call__(self, record):
        if re.search(self.pattern, record.CHROM):
            print(re.search(self.pattern, record.CHROM).start(), record.CHROM)
            return record

    def filter_name(self):
        return "%s-/%s/" % (self.name, self.pattern)

    def __repr__(self):
        return self.filter_name()
class FilterAltSubstX(VCFFilterBase):
    '''Remove variants whose only alternate allele is an X'''
    name = "badalt-X"

    @classmethod
    def customize_parser(self, parser):
        pass
    
    def __init__(self, args):
        self.threshold = args.threshold
        self.hard = args.hard

    def __call__(self, record):
        if "X" in map(str, record.ALT) and len(record.ALT) == 1 and self.hard == 0:
            return record
        elif "X" in map(str, record.ALT) and record.INFO['DP'][0] < self.threshold or len(record.ALT) == 1 and self.hard == 1:
            return record
        elif "X" in map(str, record.ALT):
            return record
        

    def filter_name(self):
        return "%s-%d-%r" % (self.name, self.threshold, self.hard)

    def __repr__(self):
        return self.filter_name()
class FilterByAltCallDepth(VCFFilterBase):
    '''Filter a VCF File by limiting variants to only those with at least X% Alt calls'''
    
    name = 'alt-call-depth'
    
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--alt-depth', type=float, default=0.4, 
            help="The minimum percentage of all calls for a locus that must be an alternative allele")
    
    def __init__(self, args):
        self.alt_depth = args.alt_depth

    def filter_name(self):
        return "%s-%.2f" % (self.name, self.alt_depth)

    def __call__(self, record):
        total_reads = sum(record.INFO["DP4"])
        ref_reads = sum(record.INFO["DP4"][:2])
        alt_reads = sum(record.INFO["DP4"][2:])

        if alt_reads / float(total_reads) >= self.alt_depth:
            return record
class FilterByReadDepth(VCFFilterBase):
    name = 'call-depth'
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--min-depth', type=int, default=5, 
            help="The minimum number of reads that must map to a location to trust a given variant call [default:5]")


    def __init__(self, args):
        self.min_depth = args.min_depth

    def filter_name(self):
        return "%s-%d" % (self.name, self.min_depth)

    def __call__(self, record):
        if record.INFO['DP'] >= self.min_depth:
            return record

#END            