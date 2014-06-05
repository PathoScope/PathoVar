import re
import os

from collections import defaultdict

## Third Party Library Imports
import vcf 
from vcf.parser import _Info

## Internal Imports
from pathovar.web.annotation_manager import EntrezAnnotationManager

from pathovar import utils
from pathovar.utils import defline_parser
from pathovar.utils import vcf_utils

def init_quality_filters(filter_args):
    filters = [filt(filter_args) for filt in vcf_utils.EXPOSED_FILTERS]
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
        self.gene_list = defaultdict(list)

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
    # Locates uncovered regions within genic or intergenic reigons
    def annotate_uncovered_genome_regions(self, genome_region_dict):
        for genome_gi, regions in genome_region_dict.items():
            opt_args = dict(verbose=self.verbose)
            for region in regions:
                genome_annotations = self.annotation_manager.get_genome_annotation(genome_gi)
                #Mutate the UncoveredRegion object
                region.location = genome_annotations.locate_snp_site(region).gid
        return genome_region_dict

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

    def get_variant_genes(self):
        genes = {}
        variant_by_gene = defaultdict(list)
        intergenic_variants = defaultdict(list)
        for var in self.annotated_variants:
            for anno in var.annotations:
                if anno.is_intergenic:
                    idents = defline_parser(var.CHROM)
                    org_name = idents.get('org_name', var.CHROM).replace('_', ' ')
                    gid = idents['gene_id']
                    intergenic_variants[org_name + ":" + gid].append(var)
                else:
                    gene_dict = {
                        "gi": anno.gid, "title": anno.title,
                        "acc": anno.accession, "strand": anno.strand,
                    }
                    genes[gene_dict['gi']] = gene_dict
                    variant_by_gene[gene_dict['gi']].append(var)
        return (genes, variant_by_gene, intergenic_variants)




class ReferenceUnparsedException(Exception):
    pass