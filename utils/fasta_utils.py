import re
from pathovar.utils import defline_parser


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
            match = re.search(r'^>(?P<defline>.+)\n', line)
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
        record = SequenceRecord(defline, sequence, self.defline_parse_func)
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
    def __init__(self, file_path, defline_parse_func = defline_parser, **opts):
        FastaParser.__init__(self, file_path, **opts)
        self.newline_between = opts.get('newline_between', False)

    def process_record(self, defline, sequence, qual):
        record = SequenceRecord(defline, sequence, self.defline_parse_func)
        record.attributes['quality'] = qual
        self.sequences.append(record)

    def parse_file(self):
        defline = ''
        sequence = ''
        qual = ''
        state = 'newline_between' if self.newline_between else 'sequence'
        for line in open(self.file_path, 'r'):
            match = re.search(r'^@(?P<defline>[^@]+)\n', line)
            if ( match and state == 'newline_between' and self.newline_between) or (match and not self.newline_between):
                if defline != '':
                    self.process_record(defline, sequence, qual)
                defline = match.groupdict()['defline']
                sequence = ''
                qual = ''
                state = 'sequence'
            elif re.search(r'\+\n', line):
                state = 'qual'
            elif line == '\n' and self.newline_between:
                state = 'newline_between'
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

    def find_uncovered_regions(self):
        undef_region = []
        last_start = None
        for ind, nucl in enumerate(self.sequence.lower()):
            if nucl == 'n':
                if last_start == None:
                    last_start = ind
            else:
                if last_start is not None:
                    undef_region.append((last_start, ind))
                    last_start = None
        return undef_region



    def __repr__(self):
        return "SequenceRecord(" + self.defline + ")"
