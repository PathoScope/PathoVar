import argparse
import glob
import os
import shutil
import requests

from pathovar.utils import config

setup_argparser = argparse.ArgumentParser()
setup_argparser.add_argument('--remove', action="store_true", default=False, help = "Remove database files")
setup_argparser.add_argument('--update', action="store_true", default=False, help = "Update database files")
setup_argparser.add_argument('--config-file', action='store', default=None, help = "Path to the configuration file to use")

def get_args():
    args = setup_argparser.parse_args()
    return args.__dict__

class SetupManager(object):
    def __init__(self, database_key, config_file = None, **opts):
        config_data = config.get_config(config_file)
        self.config_data = config_data
        self.database_data = config_data['external_databases'][database_key]
        self.database_key = database_key
        self.verbose = opts.get('verbose', False)
        self.storage_path = config_data["database_storage_directory"] + os.sep + self.database_data['storage_path']
        self.opts = opts
        self.errors = 0

        print("Storage Path: " + self.storage_path)

    def _pre_setup_dirs(self, *args, **kwargs):
        pass
    def _post_setup_dirs(self, *args, **kwargs):
        pass

    def setup_dirs(self, *args, **kwargs):
        if self.verbose: print("Setting up directories")
        try:
            os.makedirs(self.storage_path)
        except OSError, e:
            if e.errno == 17: 
                pass
            else:
                if self.verbose: print("Ignoring: ", e, e.errno)
                self.errors += 1
        self._pre_setup_dirs(*args, **kwargs)
        for key, value in self.database_data['data_urls'].items():
            if type(value) == list:
                try:
                    if self.verbose: print("Setting up %s" % os.path.join(self.storage_path, key))
                    os.makedirs(os.path.join(self.storage_path, key))
                except OSError, e:
                    if e.errno == 17: 
                        pass
                    else:
                        if self.verbose: print("Ignoring: ", e, e.errno)
                        self.errors += 1
        self._post_setup_dirs(*args, **kwargs)

    def _pre_remove_dirs(self, *args, **kwargs):
        pass
    def _post_remove_dirs(self, *args, **kwargs):
        pass

    def remove_dirs(self, *args, **kwargs):
        if self.verbose: print("Setting up directories")
        try:
            os.makedirs(self.storage_path)
        except OSError, e:
            if e.errno == 17: 
                pass
            else:
                if self.verbose: print("Ignoring: ", e, e.errno)
                self.errors += 1
        self._pre_remove_dirs(*args, **kwargs)
        for key, value in self.database_data['data_urls'].items():
            if type(value) == list:
                try:
                    if self.verbose: print("Removing %s" % os.path.join(self.storage_path, key))
                    map(os.remove, glob.glob(os.path.join(self.storage_path, key, '*')))
                    os.rmdir(os.path.join(self.storage_path, key))
                except OSError, e:
                    if e.errno == 17: 
                        pass
                    else:
                        if self.verbose: print("Ignoring: ", e, e.errno)
                        self.errors += 1
        self._post_remove_dirs(*args, **kwargs)

    def _pre_download_file(self, *args, **kwargs):
        pass

    def _post_download_file(self, *args, **kwargs):
        pass

    def _post_setup_actions(self, *args, **kwargs):
        pass

    def download_files(self, *args, **kwargs):
        for key, value in self.database_data['data_urls'].items():
            if type(value) == list:
                destination = os.path.join(self.storage_path, key)
                for url in value:
                    self._pre_download_file(url, *args, **kwargs)
                    result = self.get_file_by_url(url, destination)
                    self._post_download_file(result, *args, **kwargs)
            else: 
                self._pre_download_file(value, *args, **kwargs)
                result = self.get_file_by_url(value, self.storage_path)
                self._post_download_file(result, *args, **kwargs)

    def get_file_by_url(self, url, destination_dir):
        if self.verbose: print("Getting %s" % url)
        file_data = requests.get(url)
        try:
            file_data.raise_for_status()
        except requests.exceptions.HTTPError, e:
            print("Failed to locate resource, %r" % e)
            self.errors += 1
            return
        file_name = url.split('/')[-1]
        file_name = file_name.split('?')[0]
        with open(destination_dir + os.sep + file_name, 'wb') as data_file:
            data_file.write(file_data.content)
        unzip = os.path.splitext(file_name)
        if unzip[1] == '.zip':
            os.system('unzip -o ' + destination_dir + os.sep + file_name + ' -d' + destination_dir)
            os.remove(destination_dir + os.sep + file_name)
            file_name = file_name[:-4]
        if unzip[1] == '.gz':
            os.system('gunzip -f ' + destination_dir + os.sep + file_name)
            file_name = file_name[:-3]
        return os.path.join(destination_dir, file_name)

    def run(self, *args, **kwargs):
        if self.opts.get('remove', False):
            self.remove_dirs(*args, **kwargs)
        elif self.opts.get('update', False):
            self.remove_dirs(*args, **kwargs)
            self.setup_dirs(*args, **kwargs)
            self.download_files(*args, **kwargs)
            self._post_setup_actions(*args, **kwargs)
        else:
            self.setup_dirs(*args, **kwargs)
            self.download_files(*args, **kwargs)
            self._post_setup_actions(*args, **kwargs)
        if self.errors > 0:
            print("Errors have occurred")
            exit(1)

    def remove(self, *args, **kwargs):
         shutil.rmtree(os.path.join( self.storage_path))