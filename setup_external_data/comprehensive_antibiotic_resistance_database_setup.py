import os
import sys
import requests
import subprocess

from pathovar import get_external_databases_config, INSTALL_DIR
from pathovar.setup_external_data import SetupManager, get_args

database_data = get_external_databases_config()["comprehensive_antibiotic_resistance_database"]

class CARDSetupManager(SetupManager):
    def __init__(self, *args, **kwargs):
        SetupManager.__init__(self, database_data, "comprehensive_antibiotic_resistance_database", *args, **kwargs)


if __name__ == '__main__':
    CARDSetupManager(verbose=True, **get_args()).run()
