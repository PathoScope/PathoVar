import re
from collections import defaultdict

from pathovar.utils import defline_parser

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
        self.reference_headers = {}
        self.meta_headers = []
        self.reads = defaultdict(list)
        self.line = None

    def parse_file(self):
        with open(self.file_name, 'r') as handle:
            for line in handle:
                self.line = line
                if "@SQ" in line:
                    result = SAMHeader.parse(line)
                    self.headers.append(result)
                    self.reference_headers[result.definition] = result
                elif "@" in line[0]:
                    result = SAMHeader.parse(line)
                    self.meta_headers.append(result)
                else: 
                    result = SAMRead.parse(line)
                    self.reads[result.rname].append(result)

    def alter_reference(self, target, replacement):
        new_headers = {}
        for name, header in self.reference_headers.items():
            if re.search(target, name):
                header.definition = re.sub(target, replacement, name)
                new_headers[re.sub(target, replacement, name)] = header
            else:
                new_headers[name] = header
        new_reads = {}
        for name, reads in self.reads.items():
            if re.search(target, name):
                for read in reads:
                    read.rname = re.sub(target, replacement, name)
                new_reads[re.sub(target, replacement, name)] = reads
            else:
                new_reads[name] = reads

        self.meta_headers[0].rest = ["SO:unsorted"]
        self.reference_headers = new_headers
        self.reads = new_reads

    def remove_reference(self, target):
        new_headers = {name:header for name, header in self.reference_headers.items() if not re.search(target, name)}
        new_reads = {name:reads for name, reads in self.reads.items() if not re.search(target,name)}
        self.meta_headers[0].rest = ["SO:unsorted"]
        self.reference_headers = new_headers
        self.reads = new_reads

    def write_out(self):
        outfile = None
        if  "out" in self.opts:
            outfile = self.opts['out']
        else:
            outfile = self.file_name + '.filt.sam'
        with open(outfile, 'w') as handle:
            handle.write(str(self.meta_headers[0]) + "\n")
            for head_line in self.reference_headers.values():
                handle.write(str(head_line))
            for head_line in self.meta_headers[1:]:
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
        aln =  (fields[11:])
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
