import os

__all__ = ["samtools_snp_caller", "snp_caller_base", "snp_utils"]

def is_bam(file_name):
    name, ext = os.path.splitext(file_name)
    return ext.lower() == ".bam"