##@ Pathovar
#

import sys 
import os
import json

__all__ = ['snp_caller', 'web', 'snp_annotation', 'tests', 'setup']

def get_external_databases_config():
    return json.load(open(os.path.dirname(__file__) + '/external_databases.json', 'r'))

##
# Installing Required Packages
#
# Make sure pip, the Python Package Manager is installed. See http://www.pip-installer.org/en/latest/installing.html
#
# python -m pip install -r /path/to/pathovar/requirements.txt --user 
#