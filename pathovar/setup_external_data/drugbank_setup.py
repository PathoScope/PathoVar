from pathovar.setup_external_data import SetupManager, get_args
from pathovar.snp_annotation import drugbank_annotator

class DrugBankSetupManager(SetupManager):
    def __init__(self, *args, **kwargs):
        SetupManager.__init__(self, 'drugbank', *args, **kwargs)

    def _post_setup_actions(self, *args, **kwargs):
        blast_bin_dir = self.config_data["tool_paths"]["blast"]
        db_blast_driver = drugbank_annotator.DrugBankProteinBlastAnnotator(self.storage_path, blast_bin_dir, verbose = True)
        for driver in db_blast_driver.blast_drivers:
            if not driver.is_built():
                raise Exception("Database not built", driver)
        print("Blast Databases Created")

if __name__ == '__main__':
    DrugBankSetupManager(verbose=True, **get_args()).run()

