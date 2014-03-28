from pathovar import get_external_databases_config, INSTALL_DIR
from pathovar.setup import SetupManager, get_args

database_data = get_external_databases_config()['immune_epitope_database']

class ImmuneEpitopedatabaseSetupManager(SetupManager):
    def __init__(self, *args, **kwargs):
        SetupManager.__init__(self, database_data, 'immune_epitope_database',*args, **kwargs)

if __name__ == '__main__':
	ImmuneEpitopedatabaseSetupManager(verbose = True, **get_args()).run()