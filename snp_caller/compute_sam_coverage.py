import gzip # Maybe make the persistence file smaller?
import json
import os
import sys

from collections import defaultdict, namedtuple

from pathovar.utils import sam_utils, defline_parser

try:
    #py3k makes xrange the default range function
    range = xrange
except Exception:
    pass

class UncoveredRegion(object):
    # No need for a forking __init__ method. Instead
    # just initialize from using UncoveredRegion.from_json(**json_dict)
    @classmethod
    def from_json(cls, start, end, location=None):
        return UncoveredRegion(start, end, location)
    
    def __init__(self, start, end, location=None):
        self.start = start
        self.end = end
        self.location = location

    def to_json_safe_dict(self):
        return self.__dict__

    def __repr__(self):
        return "UncoveredRegion(start:%(start)s end:%(end)s location:%(location)s)".format(**self.__dict__)
##
# Computes the per-base coverage for a region of the genome
def extract_mean_coverage(coverage_dict, start, end):
    total_coverage = 0
    for pos in range(start, end):
        pos = str(pos)
        if pos in coverage_dict:
            total_coverage += coverage_dict[pos]
    mean_coverage = total_coverage/(end-start)
    return mean_coverage

def find_uncovered_regions(coverage_dict):
    uncovered_spans = dict()
    for genome in coverage_dict:
        last_ind = 0
        regions = []
        for pos in sorted(coverage_dict[genome]):
            if last_ind < (pos - 1):
                regions.append(UncoveredRegion(last_ind + 1, pos - 1,))
            last_ind = pos
        uncovered_spans[genome] = regions
    return uncovered_spans

# Only compatible with .sam files. Parser does not 
# work on .bam files. 
def compute_coverage(sam_file = None, sam_parser = None):
    if sam_parser is None and sam_file is not None:
        sam_parser = sam_utils.SamParser(sam_file)

    coverage_dict = defaultdict(lambda : defaultdict(int))
    uncovered_regions_dict = dict()

    for reference, reads in sam_parser.reads.items():
        gid = defline_parser(reference)['gene_id']
        genome_len = sam_parser.reference_headers[reference].fields['LN']
        coverage_dict[gid][genome_len] = 0
        for read in reads:
            for i in range(int(read.pos), int(read.pos) + len(read.seq)):
                coverage_dict[gid][i]+=1

    uncovered_regions_dict = find_uncovered_regions(coverage_dict)
    return coverage_dict, uncovered_regions_dict

def to_json_safe_dict(coverage_dict, uncovered_regions_dict):
    uncovered_regions_dict = {genome: [reg.to_json_safe_dict() for reg in regions] for
                                genome, regions in uncovered_regions_dict.items()}
    result = {"uncovered_regions_dict": uncovered_regions_dict, "coverage_dict": coverage_dict}
    return result

def from_json(json_path):
    json_dict = json.load(open(json_path))
    uncovered_regions_dict = json_dict["uncovered_regions_dict"]
    uncovered_regions_dict = {genome: [UncoveredRegion.from_json(reg['start'], reg['end'], reg['location']) for reg in regions] for
                                genome, regions in uncovered_regions_dict.items()}
    json_dict["uncovered_regions_dict"] = uncovered_regions_dict
    return json_dict

def main(sam_file = None, sam_parser = None):
    coverage_dict, uncovered_regions_dict = compute_coverage(sam_file, sam_parser)
    file_name = sam_parser.file_name if sam_parser is not None else sam_file
    result = to_json_safe_dict(coverage_dict, uncovered_regions_dict)
    outfile = os.path.splitext(file_name)[0] + '.coverage.json'
    json.dump(result, open(outfile,'w'))
    return outfile

if __name__ == '__main__':
    main(sam_file = sys.argv[1])


