import os
import sys
import requests
import subprocess

from pathovar import get_external_databases_config, INSTALL_DIR
from pathovar.setup import SetupManager, get_args

database_data = get_external_databases_config()['drugbank']

class DrugBankSetupManager(SetupManager):
    def __init__(self, *args, **kwargs):
        SetupManager.__init__(self, database_data, *args, **kwargs)

if __name__ == '__main__':
    DrugBankSetupManager(verbose=True, **get_args()).run()

