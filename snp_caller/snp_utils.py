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
    def __init__(self, file_name, out= None, **opts):
        self.file_name = file_name
        self.opts = opts
        if out != None:
            self.opts['out'] = out
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
        self.headers[0].rest = ["SO:unsorted"]

    def remove_reference(self, target):
        header = filter(lambda x: x.definition == target, self.headers)[0]
        del self.reads[header.definition]
        del header
        self.headers[0].rest = ["SO:unsorted"]

    def write_out(self):
        outfile = None
        if  "out" in self.opts:
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
genbank_id = re.compile(r'gb\|([^\|]+)\|')
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
    genbank_id_matches = genbank_id.findall(line)
    if len(genbank_id_matches) > 0:
        result['genbank_id'] = genbank_id_matches[0]
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

class FilterByAltCallDepth(VCFFilterBase):
    '''Filter a VCF File by limiting variants to only those with at least X% Alt calls'''
    
    name = 'alt-call-depth'
    
    @classmethod
    def customize_parser(self, parser):
        parser.add_argument('--percent-alt', type=float, default=0.4, 
            help="The minimum percentage of all calls for a locus that must be an alternative allele")
    
    def __init__(self, args):
        self.percent_alt = args.percent_alt

    def filter_name(self):
        return "%s-%f" % (self.name, self.percent_alt)

    def __call__(self, record):
        total_reads = sum(record.INFO["DP4"])
        ref_reads = sum(record.INFO["DP4"][:2])
        alt_reads = sum(record.INFO["DP4"][2:])

        if alt_reads / total_reads >= self.percent_alt:
            return record

## Namespace
# For compatibility with classes expecting 
# input from argparse
class Namespace:
    def __init__(self):
        pass

    def __repr__(self):
        return str(self.__dict__)

## 
# Parse a Fasta format file and provide basic filtering
# utilities to keep sequences that meet a certain criteria
class FastaParser(object):
    def __init__(self, file_path, defline_parse_func = defline_parser, **opts):
        self.file_path = file_path
        self.defline_parse_func = defline_parse_func
        self.opts = opts
        self.outfile_path = opts.get("out", file_path + '.filtered')
        self.sequences = []
        self.parsed = False

    def __iter__(self):
        if not self.parsed: self.parse_file()
        for sequence in self.sequences:
            yield sequence

    ## parse_file
    # Extracts all sequences from the fasta file. on disk, 
    # converting them into `SequenceRecords` using `process_record`
    def parse_file(self):
        defline = ''
        sequence = ''
        for line in open(self.file_path, 'r'):
            match = re.search(r'^>(?P<defline>.+\n)', line)
            if match:
                if defline != '':
                    self.process_record(defline, sequence)
                defline = match.groupdict()['defline']
                sequence = ''
            else:
                sequence += line
        self.process_record(defline, sequence)

        self.parsed = True


    ## process_record
    # Combine a defline and a sequence into a 
    # `SequenceRecord` object
    # @param defline Fasta-format defline string
    # @param sequence Single-character-to-residue string
    def process_record(self, defline, sequence):
        record = SequenceRecord(defline, sequence)
        self.sequences.append(record)

    ## filter_by_org_name
    # Filter the read sequences, keeping only those whose
    # `org_name` field matches the regular expression provided
    # @param org_name_regex Regular expresson matching an organism name
    # @sideeffect Modifies self.outfile_path to include the sanitized regular
    # expression
    def filter_by_org_name(self, org_name_regex):
        keepers = [record for record in self.sequences if re.search(org_name_regex, record.org_name)]
        self.sequences = keepers
        self.outfile_path += '.org_' + re.sub(r'[/\\:*?"<>|{}]', '_', org_name_regex)

    ## 
    # @param tax_ids_regex Regular expression to match taxonomy ids 
    def filter_by_tax_ids(self, tax_ids_regex):
        keepers = [record for record in self.sequences if re.search(tax_ids_regex, record.tax_id)]
        self.sequences = keepers
        self.outfile_path += '.tis_' + re.sub(r'[/\\:*?"<>|{}]', '_', tax_ids_regex)
    
    ## filter_by_gene_ids
    # Filter the read sequences, keeping only those whose
    # `gene_id` field is in the set of gene_ids provided
    def filter_by_gene_ids(self, gene_ids_regex):
        keepers = [record for record in self.sequences if re.search(gene_ids_regex, record.gene_id)]
        self.sequences = keepers
        self.outfile_path += '.gis_' + re.sub(r'[/\\:*?"<>|{}]', '_', gene_ids_regex)

    def filter_by_defline(self, defline_regex):
        keepers = [record for record in self.sequences if re.search(defline_regex, record.defline)]
        self.sequences = keepers
        self.outfile_path += '.defline_' + re.sub(r'[/\\:*?"<>|{}]', '_', defline_regex)

    ## 
    # Writes the remaining sequences to file in Fasta Format
    def write_output(self):
        outfile = open(self.outfile_path + '.fa', 'w')
        for record in self.sequences:
            outfile.write(record.to_fasta_format())
        outfile.close()
        return outfile.name

class FastQParser(FastaParser):
    def __init__(self, file_path, **opts):
        super(FastaParser, self).__init__(self, file_path, **opts)

    def process_record(self, defline, sequence, qual):
        record = SequenceRecord(defline, sequence)
        record.attributes['quality'] = qual
        self.sequences.append(record)

    def parse_file(self):
        defline = ''
        sequence = ''
        qual = ''
        state = 'sequence'
        for line in open(self.file_path, 'r'):
            match = re.search(r'^@(?P<defline>[^@]+\n)', line)
            if match:
                if defline != '':
                    self.process_record(defline, sequence, qual)
                defline = match.groupdict()['defline']
                sequence = ''
                qual = ''
                state = 'sequence'

            elif re.search(r'\+\n', line):
                state = 'qual'
            else:
                if state == 'sequence':
                    sequence += line
                else: 
                    qual += line
        self.process_record(defline, sequence, qual)
        self.parsed = True

    def write_output(self):
        outfile = open(self.outfile_path + '.fq', 'w')
        for record in self.sequences:
            outfile.write(record.to_fastq_format())
        outfile.close()
        return outfile.name

class SequenceRecord(object):
    def __init__(self, defline, sequence, defline_parser_func):
        self.defline = defline
        defline_fields = defline_parser_func(defline)
        self.org_name = defline_fields.get('org_name','-')
        self.tax_id = defline_fields.get('tax_id','-')
        self.gene_id = defline_fields.get('gene_id', '-')
        self.sequence = sequence
        self.attributes = defline_fields

    def to_fasta_format(self):
        entry = ">" + self.defline
        entry += self.sequence
        return entry

    def to_fastq_format(self):
        entry = "@" + self.defline
        entry += self.sequence
        entry += '+\n'
        qual = self.attributes.get('quality', None)
        if qual == None:
            qual = "!" * len(self.sequence)
        entry += qual
        return entry

    def __repr__(self):
        return "SequenceRecord(" + self.defline + ")"



