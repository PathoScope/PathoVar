##@samtools_snp_caller.py
# 
import os
import sys
import subprocess

from pathovar.snp_caller import is_bam
from pathovar.snp_caller.snp_caller_base import SNPCallerBase, SNPCallerException

## 
#
# Defines the logic for using Samtools as the SNP Caller. 
# Assumes that all of samtools, vcfutils and bcftools are 
# stored in the same directory
class SamtoolsSNPCaller(SNPCallerBase):
##Constructor
# @param opts A defaultdict(bool) that carries configuration options
# @param bin_dir The location of the `samtools` executable files
    def __init__(self, bin_dir = '', **opts):
        SNPCallerBase.__init__(self, bin_dir, **opts)
        self.intermediary_files = []

## 
#
    def get_reference_genome(self, reference_fasta, **kwargs):
        return SNPCallerBase.get_reference_genome(self, reference_fasta, **kwargs)

    def drop_missing_references_from_alignment(self, sam_file, reference_fasta, **kwargs):
        return SNPCallerBase.drop_missing_references_from_alignment(self, sam_file, reference_fasta, **kwargs)

## 
# @param self The object reference
# @param output_sam_file The path to the target .sam file
# @param source The path to the reference genome of all of the organisms of interest
    def call_snps(self, sam_file, source,  **kwargs):
        genome_path = self.get_reference_genome(source, **kwargs)
        indexed_reference_genome = self._index_reference_genome(genome_path)
        # Input may be in BAM format
        bam_file = [sam_file]
        if not is_bam(sam_file):
            if self.verbose: print("Input is .sam, filtering")
            cleaned_sam_file = self.drop_missing_references_from_alignment(sam_file, genome_path)
            bam_file = self._sam_to_bam(genome_path, *[cleaned_sam_file])
        else:
            if self.verbose: print("Input is .bam, not filtering")
        sorted_bam_files = self._sort_bam(*bam_file)
        indexed_bam_files = self._index_bam(*sorted_bam_files)
        vcf_files = self._mpileup_bam_files(genome_path, *sorted_bam_files)
        if self.opts['clean']:
            self._clean_up_files(*self.intermediary_files)
        return vcf_files[0]
## 
# Runs Samtools's faidx command on a given reference genome file in fasta format
# @param self The object reference
# @param genome_path The path to the reference genome file in fasta format
# @return The path to the indexed genome
    def _index_reference_genome(self, genome_path):
        result = os.system('%ssamtools faidx %s' % (self.bin_dir, genome_path))
        if result != 0: 
            raise SNPCallerException("An error occurred when trying to index %s" % genome_path)
        indexed_genome = genome_path + '.fai'
        return indexed_genome
## 
# Converts a an arbitrary number of .sam files to .bam files. Requires an indexed reference genome.
# @param self The object reference
# @param genome_path The path to the indexed genome fasta file
# @param *sam_files Arbitrary list of string paths to .sam files to be converted under the given indexed genome
# @return A list of paths to .bam files created from the input list of .sam files
    def _sam_to_bam(self, genome_path, *sam_files):
        bam_files = []
        for sam_file in sam_files:
            bam_file = sam_file[:-3] + 'bam'
            result = os.system('%ssamtools view -hb -t %s -S %s > %s ' % (self.bin_dir, genome_path, sam_file, bam_file))
            if result != 0:
                raise SNPCallerException("An error occurred when trying to convert %s to a BAM file" % sam_file)
            bam_files.append(bam_file)
        self.intermediary_files.extend(bam_files)
        return bam_files
## 
# Calls samtools sort over an arbitrary list of .bam files
# @param bam_files An arbitrary list of paths to .bam files to be sorted.
# @return A list of paths to the sorted .bam files.
    def _sort_bam(self, *bam_files):
        sorted_bam = []
        for bam_file in bam_files:
            sort = bam_file[:-3] + 'sorted'
            print("Sorting " + bam_file)
            result = os.system('%ssamtools sort %s %s' % (self.bin_dir, bam_file, sort))
            if result != 0:
                raise SNPCallerException("An error occurred when trying to sort %s" % bam_file)
            sorted_bam.append(sort + '.bam')
        self.intermediary_files.extend(sorted_bam)
        return sorted_bam
## 
# Calls samtools index over an arbitrary list of sorted .bam files
# @param bam_files An arbitrary list of sorted .bam files to be indexed
# @return A list of paths to the indexed .bam files
    def _index_bam(self, *bam_files):
        indexed_bam_files = []
        for bam_file in bam_files:
            print("Indexing " + bam_file)
            indexed = bam_file + '.bai'
            result = os.system('%ssamtools index %s' % (self.bin_dir, bam_file) )
            if result != 0:
                raise SNPCallerException("An error occured when trying to index %s" % bam_file)
            indexed_bam_files.append(indexed)
        self.intermediary_files.extend(indexed_bam_files)
        return indexed_bam_files
## 
# Perform SNP Calling using `samtools mpileup | bcftools view | vcfutils` 
# workflow piping. 
# @param genome_path The path to a genome that has been indexed.
# @param indexed_bam_files An arbitrary list of indexed BAM files to call SNPs on
# @param opts Arbitrary keyword dictionary to extract any options from
# @return A list of paths to the resulting .vcf files
    def _mpileup_bam_files(self, genome_path, *indexed_bam_files, **opts):
        vcf_files = []
        opts['bin_dir'] = self.bin_dir
        opts['ref_genome'] = genome_path
        if 'max_read_depth' not in opts: opts['max_read_depth'] = 100
        for bam_file in indexed_bam_files:
            opts['bam_file'] = bam_file
            opts['final_vcf'] = bam_file[:-11] + '.vcf'
            opts['intermediary'] = bam_file + '.raw.bcf'
            opts['mpileup-m'] = 3
            opts['mpileup-F'] = 0.0002
            opts['consensus-fq'] = ".cns.fq"
            
            if self.verbose: print('{bin_dir}samtools mpileup -uD -m {mpileup-m} -F {mpileup-F} -f {ref_genome} {bam_file} | {bin_dir}bcftools view -bvcg - > {intermediary}'.format(**opts))

            self.snp_call_process = subprocess.Popen('{bin_dir}samtools mpileup -uD -m {mpileup-m} -F {mpileup-F} -f {ref_genome} {bam_file} | {bin_dir}bcftools view -bvcg - > {intermediary}'.format(**opts), 
                stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
            
            # Let this process run independently in the background
            consensus_process = subprocess.Popen('{bin_dir}samtools mpileup -uf {ref_genome} {bam_file} | {bin_dir}bcftools view -cg - | vcfutils.pl vcf2fq > {final_vcf}{consensus-fq}'.format(**opts),
                stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
            
            result = self.snp_call_process.wait()

            #result = os.system('{bin_dir}samtools mpileup -uD -m {mpileup-m} -F {mpileup-F} -f {ref_genome} {bam_file} | {bin_dir}bcftools view -bvcg - > {intermediary}'.format(**opts))
            if result != 0:
                raise SNPCallerException("An error occurred during samtools mpileup | bcftools view for %s" % opts['bam_params'])

            if self.verbose: print('{bin_dir}bcftools view {intermediary} | {bin_dir}vcfutils.pl varFilter -D{max_read_depth} > {final_vcf}'.format(**opts))
            self.snp_call_process = subprocess.Popen('{bin_dir}bcftools view {intermediary} | {bin_dir}vcfutils.pl varFilter -D{max_read_depth} > {final_vcf}'.format(**opts),
                stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)

            result = self.snp_call_process.wait()
            result = os.system('{bin_dir}bcftools view {intermediary} | {bin_dir}vcfutils.pl varFilter -D{max_read_depth} > {final_vcf}'.format(**opts))
            if result != 0:
                raise SNPCallerException("An error occurred during bcftools view | vcfutils.pl varFilter for %s" % bam_file)
            
            if self.verbose: print('{bin_dir}samtools mpileup -uf {ref_genome} {bam_file} | {bin_dir}bcftools view -cg - | vcfutils.pl vcf2fq > {final_vcf}{consensus-fq}'.format(**opts))
            #result = consensus_process.wait()
            #result = os.system('{bin_dir}samtools mpileup -uf {ref_genome} {bam_file} | {bin_dir}bcftools view -cg - | vcfutils.pl vcf2fq > {final_vcf}{consensus-fq}'.format(**opts))
            #if result != 0:
            #    raise SNPCallerException("An error occurred during vcfutils.pl vcf2fq %s" % bam_file)

            vcf_files.append(opts["final_vcf"])
            self.intermediary_files.append(opts['intermediary'])
        return vcf_files

## 
# A simple helper function for removing arbitrary lists of intermediary 
# files and indices that will be called when --clean is used
# @param files_to_remove Arbitrary list of file paths to be removed
    def _clean_up_files(self, *files_to_remove):
        for file_path in files_to_remove:
            print("Removing %s" % file_path)
            os.remove(file_path)
