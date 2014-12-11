import os
import subprocess
import argparse

from pathovar.utils import config

def main():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-a","--all", action='store_true', required=False, help="Install/Update all databases")
    arg_parser.add_argument("-u", "--update", action="store_true", required=False, help="Update installed databases")
    arg_parser.add_argument("-r", "--remove", action="store_true", default = False, required=False, help="Remove the datbases, don't install or update them")
    arg_parser.add_argument("-c", "--config", action="store", default = None, help = "The configuration file to use. If no file exists at that path, \
        it will be created. If not given, the default configuration is used")
    arg_parser.add_argument("databases", nargs="*", help="The names of databases to operate on")
    args = arg_parser.parse_args()
    print(args)
    databases = set(args.databases)

    if args.config is None:
        args.config = config.DEFAULT_CONFIG_PATH

    database_config_data = config.get_config(args.config)

    if("version" not in database_config_data or database_config_data['version'] != config.DEFAULT_CONFIG["version"]):
        print("Backing up old configuration file.")
        os.rename(args.config, args.config + ".old")
        print("Writing new configuration file.")
        config.init_config(args.config)
        print("Please update %s with your custom settings.\nOld settings are in %s" % (args.config, args.config + ".old"))

        databases = set([name for name, values in database_config_data['external_databases'].items() if name != "version" and values['enabled']])
        database_config_data = database_config_data = config.get_config(args.config)

    if args.all or args.databases == []:
        databases = set([name for name, values in database_config_data['external_databases'].items() if name != "version"])

    for database in databases:
        print("Handling %s --------------------------------------" % database)
        try:
            # Get the external database data and access it's setup python script.
            config_data = database_config_data["external_databases"][database]
            if config_data['setup_script'] != None:
                cmd = 'python %s' % config_data['setup_script']
                
                # Pass along action selection
                if(args.update):
                    cmd += ' --update'
                elif(args.remove):
                    cmd += ' --remove'

                # Link to the given configuration file
                cmd += ' --config-file ' + args.config

                result = subprocess.Popen(cmd, shell=True).wait()
                if(result != 0):
                    print("An error prevented this database from being set up automatically.")
                config_data['enabled'] = (result == 0)
                database_config_data[database] = config_data
            else:
                print("No setup script found for %s" % database) 
        except KeyError:
            print("No entry for %s" % database)

    config.write_config(args.config, database_config_data)

if __name__ == '__main__':
    main()