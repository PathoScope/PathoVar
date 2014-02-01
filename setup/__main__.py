import json
import sys
import subprocess

import pathovar

database_config_data = pathovar.get_external_databases_config()

databases = []

if len(sys.argv) > 1:
    databases = sys.argv[1:]

else:
    databases = database_config_data


for database in databases:
    print("Handling %s" % database)
    try:
        config_data = database_config_data[database]
        if config_data['setup_script'] != None:
            subprocess.call('python %s %s %s' % (config_data['setup_script'], config_data['storage_path'], " ".join(config_data["data_urls"])), shell=True)
            config_data['enabled'] = True
            database_config_data[database] = config_data
        else:
            print("No setup script found for %s" % database) 

    except KeyError, e:
        print("No entry for %s" % database)

pathovar.update_external_databases_config(database_config_data)