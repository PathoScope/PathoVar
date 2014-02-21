##@ Pathovar
#

import sys 
import os
import json

__all__ = ['snp_caller', 'web', 'snp_annotation', 'tests', 'setup']

DEFAULT_EXTERNAL_DATABASE_CONFIG = {
    "version": "0.2b",
    "comprehensive_antibiotic_resistance_database": {
        "data_urls" : {
            "nucleotide" : [
            "http://arpcard.mcmaster.ca/blast/db/nucleotide/AR-genes.fa.gz", 
            "http://arpcard.mcmaster.ca/blast/db/nucleotide/AT-genes.fa.gz", 
            "http://arpcard.mcmaster.ca/blast/db/nucleotide/ABS-genes.fa.gz"
            ], 
            "protein":[
            "http://arpcard.mcmaster.ca/blast/db/protein/AR-polypeptides.fa.gz",
            "http://arpcard.mcmaster.ca/blast/db/protein/AT-polypeptides.fa.gz", 
            "http://arpcard.mcmaster.ca/blast/db/protein/ABS-polypeptides.fa.gz"
            ],
            "other":[
                "http://arpcard.mcmaster.ca/obo-download/aro.obo"
            ]},
        "storage_path": "databases/comprehensive_antibiotic_resistance_database/",
        "setup_script": "setup/comprehensive_antibiotic_resistance_database_setup.py",
        "appropriate_organisms": [
            "*"
        ],
        "enabled": False
    },
    "patric": {
        "data_urls": [
            "ftp://ftp.patricbrc.org/patric2/genomes/"
        ],
        "storage_path": "databases/patric/",
        "setup_script": None,
        "appropriate_organisms": [
            "bacteria"
        ],
        "enabled": False
    },
    "drugbank": {
        "data_urls": {
            "nucleotide": [
                "http://www.drugbank.ca/system/downloads/current/sequences/gene/all_target.fasta.zip",
                "http://www.drugbank.ca/system/downloads/current/sequences/gene/all_enzyme.fasta.zip",
                "http://www.drugbank.ca/system/downloads/current/sequences/gene/all_transporter.fasta.zip",
                "http://www.drugbank.ca/system/downloads/current/sequences/gene/all_carrier.fasta.zip"
            ],
            "protein":[
                "http://www.drugbank.ca/system/downloads/current/sequences/protein/all_target.fasta.zip",
                "http://www.drugbank.ca/system/downloads/current/sequences/protein/all_enzyme.fasta.zip",
                "http://www.drugbank.ca/system/downloads/current/sequences/protein/all_transporter.fasta.zip",
                "http://www.drugbank.ca/system/downloads/current/sequences/protein/all_carrier.fasta.zip"
            ]
        },
        "storage_path": "databases/drugbank/",
        "setup_script": "setup/drugbank_setup.py",
        "appropriate_organisms": [
            "*"
        ],
        "enabled": False
    }, 
    "immune_epitope_database": {
         "data_urls":{
            "cell_assays": [
                "http://www.iedb.org/doc/tcell_compact.zip",
                "http://www.iedb.org/doc/bcell_compact.zip"            
            ]

         },
         "storage_path": "databases/immune_epitope_database/",
         "setup_script": "setup/immune_epitope_database_setup.py",
         "appropriate_organisms": [
             "*"
         ],
         "enabled": False
    }
}

def get_external_databases_config():
    global DEFAULT_EXTERNAL_DATABASE_CONFIG
    try:
        conf = json.load(open(os.path.dirname(__file__) + '/external_databases.json', 'r'))
    except Exception, e:
        print(e)
        print("Writing Fresh external_databases.json config")
        update_external_databases_config(DEFAULT_EXTERNAL_DATABASE_CONFIG)
        conf = get_external_databases_config()
    if "version" not in conf or conf['version']:
        print("Your external database configuration is outdated, please rerun `python -m pathovar.setup` to update it!")
    return conf

def update_external_databases_config(config_dict):
    json.dump(config_dict, open(os.path.dirname(__file__) + '/external_databases.json', 'w'))

INSTALL_DIR = os.path.dirname(__file__)




##
# Installing Required Packages
#
# Make sure pip, the Python Package Manager is installed. See http://www.pip-installer.org/en/latest/installing.html
#
# python -m pip install -r /path/to/pathovar/requirements.txt --user 
#
