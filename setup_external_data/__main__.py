import json
import os
import sys
import subprocess
import argparse
import pathovar


def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-a","--all", action='store_true', required=False, help="Install/Update all databases")
    arg_parser.add_argument("-u", "--update", action="store_true", required=False, help="Update installed databases")
    arg_parser.add_argument("-r", "--remove", action="store_true", default = False, required=False, help="Remove the datbases, don't install or update them")
    arg_parser.add_argument("-c", "--config", action="store", default = pathovar.INSTALL_DIR, help = "The configuration file to use")
    arg_parser.add_argument("databases", nargs="*", help="The names of databases to operate on")
    args = arg_parser.parse_args()
    print(args)
    databases = set(args.databases)

    # Attach args.config here
    if(not os.path.exists(args.config)):
        pass
    database_config_data = pathovar.get_external_databases_config()
    
    if("version" not in database_config_data or 
        database_config_data['version'] != pathovar.DEFAULT_EXTERNAL_DATABASE_CONFIG["version"]):
        print("Updating external database configuration file...")
        pathovar.update_external_databases_config(pathovar.DEFAULT_EXTERNAL_DATABASE_CONFIG)
        databases = set([name for name, values in database_config_data.items() 
            if name != "version" and values['enabled']])
        database_config_data = pathovar.get_external_databases_config()


    if args.all or args.databases == []:
        databases = set([name for name, values in database_config_data.items() 
            if name != "version"])

    for database in databases:
        print("Handling %s" % database)
        try:
            config_data = database_config_data[database]
            if config_data['setup_script'] != None:
                cmd = 'python %s' % config_data['setup_script']
                if(args.update):
                    cmd += ' --update'
                elif(args.remove):
                    cmd += ' --remove'

                # Attach args.config here
                cmd += ' --config-file ' + args.config

                result = subprocess.check_call(cmd, shell=True)
                config_data['enabled'] = (result == 0)
                database_config_data[database] = config_data
            else:
                print("No setup script found for %s" % database) 
        except subprocess.CalledProcessError, e:
            config_data['enabled'] = False
            database_config_data[database] = config_data
        except KeyError, e:
            print("No entry for %s" % database)

    pathovar.update_external_databases_config(database_config_data)

if __name__ == '__main__':
    main()