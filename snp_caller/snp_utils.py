##@ snp_utils.py
# A collection of utilities for exploring SNP data and data formats

# Standard Library Imports
import os
import sys
import re
from collections import defaultdict

# 3rd Party Imports
import vcf
from vcf.parser import _Filter
from vcf.filters import Base as VCFFilterBase

## SpeciesExtractNameSAMParser
# 
class SpeciesExtractNameSAMParser(object):
    def __init__(self, file_name):
        self.file_name = file_name
        self.species = defaultdict(list)
        self.line_buffer = []

    def parse_file(self):
        handle = open(self.file_name, 'rb')
        lines = handle.readlines()
        for line in lines:
            #print("LINE: " + line)
            if re.match(r"^##", line):
                self.buffer.append(line)
            if re.match(r"^@", line):
                self.buffer.append(line)
            fields = line.split('\t')
            #print("FIELDS: " + str(fields))
            #raw_input("NEXT?")
            self.species[fields[2]].append(fields)
        return self.species

    def filter_species(self):
        species_data = self.species.items()
        species_data.sort(key=lambda x:len(x[1]), reverse=True)
        return species_data[:5]

class SAMParser(object):
    def __init__(self, file_name, **opts):
        self.file_name = file_name
        self.opts = opts
        self.headers = []
        self.reads = defaultdict(list)
        self.line = None

    def parse_file(self):
        with open(self.file_name, 'r') as handle:
            for line in handle:
                self.line = line
                if "@" in line:
                    result = SAMHeader.parse(line)
                    self.headers.append(result)
                else: 
                    result = SAMRead.parse(line)
                    self.reads[result.rname].append(result)

    def alter_reference(self, target, replacement):
        header = filter(lambda x: x.definition == target, self.headers)[0]
        reads = self.reads[header.definition]
        for read in reads:
            read.rname = replacement
        header.definition = replacement
        header[0].rest = ["SO:unsorted"]

    def write_out(self):
        outfile = None
        if "out" in self.opts:
            outfile = self.opts['out']
        else:
            outfile = self.file_name + '.filt.sam'
        with open(outfile, 'w') as handle:
            for head_line in self.headers:
                handle.write(str(head_line))
            for rname in self.reads:
                for read in self.reads[rname]:
                    handle.write(str(read))

class SAMHeader(object):
    @classmethod
    def parse(self, line):
        fields = line.split('\t')
        result = SAMHeader(fields)
        return result

    def __init__(self, fields):
        self.tag = fields[0]
        self.definition = fields[1]
        if len(fields) > 2:
            self.rest = fields[2:]
        else: 
            self.rest = []

    def __repr__(self):
        rest = '\t'.join(self.rest)
        rep = "{tag}\t{definition}\t{rest}".format(**{"tag": self.tag, "definition": self.definition, "rest": rest})
        return rep

class SAMRead(object):
    @classmethod
    def parse(self, line):
        fields = line.split("\t")
        aln =  (fields[12:])
        fields = fields[:11]
        fields.append(aln)
        result = SAMRead(*fields)
        return result

    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, aln):
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        self.aln = '\t'.join(aln)

    def __repr__(self):
        rep = "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\t{aln}".format( **self.__dict__ )
        return rep



##
gene_id = re.compile(r'gi\|([^\|]+)\|')
##
tax_id = re.compile(r'ti\|([^\|]+)\|')
##
org_name = re.compile(r'org[^\|]*\|([^\|]+)\|')

##
#
def defline_parser(line):
    result = {}
    org_name_matches = org_name.findall(line)
    if len(org_name_matches) > 0:
        result['org_name'] = org_name_matches[0]
    tax_id_matches = tax_id.findall(line)
    if len(tax_id_matches) > 0:
        result['tax_id'] = tax_id_matches[0]
    gene_id_matches = gene_id.findall(line)
    if len(gene_id_matches) > 0:
        result['gene_id'] = gene_id_matches[0]
    if len(result.keys()) == 0:
        raise DeflineUnparsedException("Could Not Parse Defline: %s" % line)
    return result

class DeflineUnparsedException(Exception):
    pass

##
#
def top_five_species_in_sam_file(file_path):
    species_extract_parser = SpeciesExtractNameSAMParser(file_path)
    species_extract_parser.parse_file()
    top_five_species = species_extract_parser.filter_species()
    result = [defline_parser(spec[0]) for spec in top_five_species]
    return result

def defline_mangler(line, match):
    raise NotImplementedError()

##
#
def report_snps_to_console(file_path):
    handle = open(file_path, 'r')
    reader = vcf.VCFReader(handle)
    for variant in reader:
        print(variant)

## filter_vcf
# Based on vcf_filter.py from PyVCF. Compatible with instantiated
# PyVCF _Filter objects.
# @param drop If True, keep sequences that pass a filter, else exclude sequences that pass the filter
def filter_vcf(file_path, filters, drop = True, short_circuit = False):
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
        if not drop:
            output_record = not output_record
        if output_record:
            # use PASS only if other filter names appear in the FILTER column
            if record.FILTER is None and not drop: record.FILTER = 'PASS'
            output.write_record(record)

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

## Namespace
# For compatibility with classes expecting 
# input from argparse
class Namespace:
    def __init__(self):
        pass

    def __repr__(self):
        return str(self.__dict__)