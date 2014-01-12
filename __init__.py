import sys 
import os
import json

__all__ = ['snp_caller', 'web', 'snp_annotation', 'tests', 'setup']

class PathoVar(object):
    def __init__(self, params):
        pass

def get_external_databases_config():
    return json.load(open(os.path.dirname(__file__) + '/external_databases.json', 'r'))