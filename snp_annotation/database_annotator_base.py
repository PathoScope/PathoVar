import os
import sys
import subprocess

from pathovar.utils import fasta_utils

EXT_TO_PARSER = {
    '.fa'  : fasta_utils.FastaParser,
    '.fq'  : fasta_utils.FastQParser,
    #'.vcf' = 
}

class ExternalDatabaseBase(object):
    def __init__(self, file_path_dict, **opts):
        self.file_paths = file_path_dict
        self.opts = opts

class ExternalDatabaseFile(object):
    def __init__(self, path, parser_class = None, **opts):
        self.path = path
        self.parser_class = parser_class
        if parser_class == None:
            ext = os.path.splitext(path)[1]
            if ext in EXT_TO_PARSER:
                self.parser_class = EXT_TO_PARSER[ext]
        self.parser = None
        self.opts = opts
        self.entries = []

    def extract(self):
        if self.parser_class == None: raise Exception("Parser Class Not Found")
        self.parser = self.parser_class(self.path, **self.opts)
        self.entries = [entry for entry in self.parser]


