import copy
import json
import os

from collections import defaultdict

from pathovar import INSTALL_DIR

# For later migrating away from config only describing external databases hard-coded into the install dir
def get_config(path = None, alert = False):
    conf = dict()
    try:
        if path is None or path is False:
            if alert: print("Using default configuration file")
            conf = get_config(DEFAULT_CONFIG_PATH)
        else:
            conf = json.load(open(path))
    except Exception, e:
        if alert:
            print(e)
            print("Trying to write fresh configuration file.")
        try:
            init_config(path)
            conf = get_config(path)
        except Exception, e:
            if alert:
                print(e)
                print("Using default configuration file")
            conf = get_config(DEFAULT_CONFIG_PATH)
    if (("version" not in conf) or (conf['version'] != DEFAULT_CONFIG['version'])) and alert:
        print("Your configuration is outdated, please rerun `python -m pathovar.setup` to update it!")
    conf_defaultdict = defaultdict(lambda : None)
    conf_defaultdict.update(conf)
    return conf_defaultdict

def write_config(path = None, data = None):
    if data is None:
        data = DEFAULT_CONFIG
    defaulted = False
    try:
        if path is None or path is False:
            path = DEFAULT_CONFIG_PATH
            defaulted = True
        json.dump(data, open(path, 'w'), sort_keys = True, indent = 4)
    except OSError, e:
        print(e)
        print("Writing configuration file failed. Terminating!")
        if defaulted: print("Attempted to update install-default configuration file")
        exit(-1)

def init_config(path):
    config = copy.deepcopy(DEFAULT_CONFIG)
    dir_name = os.path.dirname(path)
    config['cache_directory'] = os.path.abspath(os.path.join(dir_name, "pathovar_annotation_cache"))
    config['database_storage_directory'] = os.path.abspath(os.path.join(dir_name, "pathovar_databases"))
    write_config(path, config)

def load_param(arg, conf, default = None):
    val = default
    if arg is not None:
        val = arg
    elif conf is not None:
        val = conf
    return val


DEFAULT_CONFIG_NAME = "pathovar.conf.json"

DEFAULT_CONFIG_PATH = os.path.join(INSTALL_DIR, DEFAULT_CONFIG_NAME)

DEFAULT_CONFIG = {
    "version": "0.4.1",
    "cache_directory": "./pathovar_annotation_cache/",
    "database_storage_directory": "./pathovar_databases/",
    "tool_paths": {
        "samtools": '',
        "snpeff": '',
        "blast": ''
    },
    "heuristic_parameters": {
        "var_score_dict": {
            "LOW": 0.1, 
            "UNKNOWN": 1,
            "MODERATE": 1, 
            "MODIFIER": 1, 
            "HIGH": 2
        },
        "blast_max": 20,
        "snp_max": 20,
        "coverage_max": 20,
        "blast_value": 2, 
        "coverage_value": 1
    },
    "filter_parameters" :{
        "alt_depth": 0.4,
        "min_depth": 5,
        "min_mq": 20,
        "min_qual": 20,
        "ref_vcfs": None,
        "intersection": False,
    },
    "external_databases":{
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
            "storage_path": "/comprehensive_antibiotic_resistance_database/",
            "setup_script": INSTALL_DIR + os.sep + "setup_external_data/comprehensive_antibiotic_resistance_database_setup.py",
            "appropriate_organisms": [
                "*"
            ],
            "enabled": False
        },
        # "patric": {
        #     "data_urls": [
        #         "ftp://ftp.patricbrc.org/patric2/genomes/"
        #     ],
        #     "storage_path": "databases/patric/",
        #     "setup_script": None,
        #     "appropriate_organisms": [
        #         "bacteria"
        #     ],
        #     "enabled": False
        # },
        "drugbank": {
        "data_urls": {
            "protein":[
                "https://sites.google.com/site/mobiuskleinscripthost/data/all_carrier.fasta.gz?attredirects=0&d=1",
                "https://sites.google.com/site/mobiuskleinscripthost/data/all_enzyme.fasta.gz?attredirects=0&d=1",
                "https://sites.google.com/site/mobiuskleinscripthost/data/all_target.fasta.gz?attredirects=0&d=1",
                "https://sites.google.com/site/mobiuskleinscripthost/data/all_transporter.fasta.gz?attredirects=0&d=1"
            ]
        },
        "storage_path": "/drugbank/",
        "setup_script": INSTALL_DIR + os.sep + "setup_external_data/drugbank_setup.py",
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
             "storage_path": "/immune_epitope_database/",
             "setup_script": INSTALL_DIR + os.sep + "setup_external_data/immune_epitope_database_setup.py",
             "appropriate_organisms": [
                 "*"
             ],
             "enabled": False
        }
    }
}



if not os.path.exists(DEFAULT_CONFIG_PATH):
    conf = copy.deepcopy(DEFAULT_CONFIG)
    conf['cache_directory'] = INSTALL_DIR + os.sep + "pathovar_annotation_cache/"
    conf['database_storage_directory'] = INSTALL_DIR + os.sep + "pathovar_databases/"
    write_config(DEFAULT_CONFIG_PATH, conf)
