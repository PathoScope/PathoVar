import re
import os

## Third Party Library Imports
import vcf 
from vcf.parser import _Info

## Internal Imports
from pathovar.web.annotation_manager import EntrezAnnotationManager

from pathovar import utils
from pathovar.utils import vcf_utils

QUAL_FILTERS = [vcf_utils.FilterByAltCallDepth, vcf_utils.FilterByReadDepth]

def init_quality_filters(filter_args):
    filters = [filt(filter_args) for filt in QUAL_FILTERS]
    return filters

class VariantLocator(object):
    def __init__(self, vcf_file, **opts):
        self.vcf_file = vcf_file
        self.opts = opts
        self.verbose = opts.get("verbose", False)

        self.reader = vcf.VCFReader(open(vcf_file, 'r'))
        self.variants = vcf_utils.filter_vcf_in_memory(self.reader, init_quality_filters(self.opts['filter_args']), keep = True)
        self.annotated_variants = []

        self.annotation_manager = opts.get("annotation_manager", EntrezAnnotationManager(**opts))

        self.reader.infos['GENE'] = _Info('GENE', 1, "String", "Gene containing this variant")

    ##
    #
    def annotate_snp(self, variant):
        ids = utils.defline_parser(variant.CHROM)
        annotations = []
        opt_args = dict(verbose=self.verbose)
        if "gene_id" in ids:
            genome_annotations = self.annotation_manager.get_genome_annotation(ids['gene_id'])
            annotations.append(genome_annotations.locate_snp_site(variant))
        else:
            raise ReferenceUnparsedException("Could not parse GID from CHROM entry, %s" % str(variant))

        anno_variant = vcf_utils.AnnotatedVariant(variant, annotations)
        return anno_variant

    ##
    #
    def annotate_all_snps(self):
        if self.verbose: print("Annotating All Variant Sites.")
        for variant in self.variants:
            anno_variant = self.annotate_snp(variant)
            self.annotated_variants.append(anno_variant)

    def write_annotated_vcf(self):
        if self.verbose: print("Writing annotated .vcf file.")
        output_file = self.vcf_file[:-4] + '.anno.vcf'
        writer = vcf.Writer(open(output_file, 'w'), self.reader)
        for variant in self.annotated_variants:
            writer.write_record(variant)
        writer.close()
        return output_file

class ReferenceUnparsedException(Exception):
    pass