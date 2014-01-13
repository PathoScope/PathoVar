import json
import sys
import subprocess

import pathovar

database_config_data = pathovar.get_external_databases_config()

for database in database_config_data:
    config_data = database_config_data[database]
    print("Handling %s" % database)
    if config_data['setup_script'] != None:
        subprocess.call('python %s %s %s' % (config_data['setup_script'], config_data['storage_path'], " ".join(config_data["data_urls"])), shell=True)
    else:
        print("No setup script found for %s" % database) 
