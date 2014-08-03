from pathovar.utils import config
from pathovar.setup_external_data import SetupManager, get_args

class ImmuneEpitopedatabaseSetupManager(SetupManager):
    def __init__(self, *args, **kwargs):
        SetupManager.__init__(self, 'immune_epitope_database',*args, **kwargs)

if __name__ == '__main__':
    ImmuneEpitopedatabaseSetupManager(verbose = True, **get_args()).run()