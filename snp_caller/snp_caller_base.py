import sys
import os
import abc

from pathovar.utils import fasta_utils
from pathovar.utils import sam_utils

class SNPCallerBase(object):
    '''Generic SNP Caller Interface class. This class does not 
directly implement any real functionality. Attempting to use any
will raise a *NotImplementedError*

    '''
    
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, bin_dir, **opts):
        '''Constructs a SNPCaller instance to make system calls 
to a third party SNP Calling program. 
@param bin_dir The location of the executables for this program
@param opts Arbitrary keyword argument configurations
'''
        self.bin_dir = bin_dir
        if self.bin_dir != '':
            self.bin_dir += os.sep
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.snp_calling_process = None
        self.consensus_sequence_process = None
        self.reference_genomes = None
        self.sam_file = None
        self.sam_parser = None


## get_reference_genome
# Retrieves and indexes the appropriate reference genome.
# @return string path to indexed reference genome file 
    @abc.abstractmethod
    def get_reference_genome(self,source, **kwargs):
        tax_ids_reg = kwargs.get("tax_ids_reg", None)
        org_names_reg = kwargs.get("org_names_reg", None) 
        gene_ids_reg = kwargs.get("gene_ids_reg", None) 
        keep_all = kwargs.get("keep_all", None) 
        if tax_ids_reg or org_names_reg or gene_ids_reg or not keep_all:
            parser = fasta_utils.FastaParser(source)
            parser.parse_file()
            if tax_ids_reg:
                parser.filter_by_tax_ids(tax_ids_reg)
            if org_names_reg:
                parser.filter_by_org_name(org_names_reg)
            if gene_ids_reg:
                parser.filter_by_gene_ids(gene_ids_reg)

            if not keep_all:
                parser.filter_by_defline(r"(complete genome)|(plasmid.*(?!encoded).*complete sequence)", ".genomes")

            filtered_source = parser.write_output()
            if os.path.getsize(filtered_source) == 0:
                raise SNPCallerException("Filtered Reference Genome File contains no genomes! Double check your filtering constraints")
            self.reference_genomes = filtered_source
            return filtered_source
        else:
            if os.path.getsize(source) == 0:
                raise SNPCallerException("Reference Genome File is empty, are you sure this is the right file? %s" % source)
            self.reference_genomes = filtered_source
            return source

    def drop_missing_references_from_alignment(self, sam_file, reference_genomes, **kwargs):
        fasta_parser = fasta_utils.FastaParser(reference_genomes) if kwargs.get('fasta_parser', None) is None \
            else kwargs.get('fasta_parser', None)
        sam_parser = sam_utils.SamParser(sam_file)
        
        sequence_gids = [seq.gid for seq in fasta_parser]
        gids_as_regex = r'gi\|' + r'|'.join(sequence_gids)
        sam_parser.remove_reference(gids_as_regex, keep = True)
        filtered_sam_file = sam_parser.write_output()
        if len(sam_parser.reference_headers) == 0:
            raise SNPCallerException("The kept reference genomes don't appear in the input Sam File. Are you sure this is the right file? %s" % sam_file)
        self.sam_file = filtered_sam_file
        self.sam_parser = sam_parser
        return filtered_sam_file


## Initiates the SNP Calling pipeline.
# @param output_sam_file A string path to the .sam file generated by PathoID
# @param tax_id integer NCBI taxonomy id
# @param org_name string organism name
# @param source string reference to reference genome location
# @param kwargs kwargs
# @return string path to .vcf or .bcf file
    @abc.abstractmethod
    def call_snps(self, output_sam_file, source, **kwargs):
        raise NotImplementedError()

class SNPCallerException(Exception):
    '''Exception class representing programmatic errors while SNP Calling'''
    pass