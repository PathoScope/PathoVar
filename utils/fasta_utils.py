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

##
# As FastaParser, but for the FastQ File format. Works for both
# 
class FastQParser(FastaParser):
    def __init__(self, file_path, defline_parse_func = defline_parser, **opts):
        FastaParser.__init__(self, file_path, **opts)
        self.multiline_mode = opts.get('multiline_mode', False)

    def process_record(self, defline, sequence, qual):
        record = SequenceRecord(defline, sequence, self.defline_parse_func)
        record.attributes['quality'] = qual
        self.sequences.append(record)

    def parse_file(self):
        defline = ''
        sequence = ''
        qual = ''
        count = 0
        state = 'start'
        for line in open(self.file_path, 'r'):
            match = re.search(r'^@(?P<defline>[^@]+)\n', line)
            if state == 'qual' and count > 0:
                qual += line
                count-=1
            elif match:
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
                    count += 1

        self.process_record(defline, sequence, qual)
        self.parsed = True

    def write_output(self):
        outfile = open(self.outfile_path + '.fq', 'w')
        for record in self.sequences:
            outfile.write(record.to_fastq_format())
        outfile.close()
        return outfile.name

class SequenceRecord(object):
    def __init__(self, defline, sequence, defline_parser_func, **attributes):
        self.defline = defline
        defline_fields = defline_parser_func(defline)
        self.org_name = defline_fields.get('org_name','-')
        self.tax_id = defline_fields.get('tax_id','-')
        self.gene_id = defline_fields.get('gene_id', '-')
        self.sequence = sequence
        self.attributes = dict(defline_fields, **attributes)

    def to_fasta_format(self):
        entry = ">" + self.defline + '\n'
        entry += self.sequence + '\n'
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

MUT_SYM = "$"
class MutatedSequenceRecord(SequenceRecord):
    def __init__(self, defline, sequence, defline_parser_func = defline_parser, **attributes):
        SequenceRecord.__init__(self, defline, sequence, defline_parser_func, **attributes)
        self.mutated_sequence = list(self.sequence)
        self.mutated_indices = []

    def sweep_mutations(self):
        for var in self.attributes['variants']:
            print("Sequence Start: ", self.attributes['start'])
            print("Variant Start: ", var["start"])
            expected_ref = str(var['ref'])
            primary_variant = str(var['alts'][0])
            print("Varinat: ", primary_variant)
            span_variant = var['end'] - var['start']
            start_offset = (var['start']-self.attributes['start'])
            # Handle variants that overlap the start of the sequence
            if start_offset < 0:
                print("Negative start_offset caught: ", start_offset)
                # Exclude the section of the reference that is outside of the gene
                expected_ref = expected_ref[(-1*start_offset):]
                primary_variant = primary_variant[(-1*start_offset):]
                span_variant += start_offset
                start_offset = 0

            # Handle variants that overlap the end of the sequence
            if start_offset + span_variant > len(self.sequence):
                print("Variant off the gene")
                span_variant -= var['end'] - self.attributes['end']
                expected_ref = expected_ref[:span_variant]
                primary_variant = primary_variant[:span_variant]

            print("start_offset: %d, span_variant: %d" % (start_offset, span_variant))
            print(''.join(self.mutated_sequence[(start_offset):(start_offset + span_variant)]), expected_ref, primary_variant, len(expected_ref) - len(primary_variant))
            assert ''.join(self.mutated_sequence[(start_offset):(start_offset + span_variant)]) == expected_ref

            # Double Check with accumulator
            insertion = []
            for index, position in enumerate(range(start_offset, start_offset + span_variant)):
                print("Pos: ", index, position)
                alternate_base = None
                if index >= len(primary_variant):
                    alternate_base = ''
                elif index == span_variant and len(primary_variant) > span_variant:
                    alternate_base = primary_variant[index:]
                else:
                    alternate_base = primary_variant[index]

                # Label each site that was mutated with a $. 
                if len(alternate_base) > 1:
                    alternate_base = MUT_SYM + MUT_SYM.join(alternate_base.split())
                else:
                    alternate_base = MUT_SYM + alternate_base
                insertion.append(alternate_base)
                self.mutated_sequence[position] = alternate_base
            print(self.mutated_sequence[start_offset:(start_offset+span_variant)], insertion)
            assert insertion == self.mutated_sequence[start_offset:(start_offset+span_variant)]
        # End of applying mutations 

        # Merge the list back into a string, now each mutation position is preceded by $
        self.mutated_sequence = ''.join(self.mutated_sequence)
        while(True):
            next_pos = self.mutated_sequence.find(MUT_SYM)
            if next_pos == -1:
                break
            self.mutated_sequence = self.mutated_sequence[:next_pos] + self.mutated_sequence[next_pos+1:]
            self.mutated_indices.append(next_pos)
        print(self.mutated_indices)

