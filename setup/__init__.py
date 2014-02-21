import os
import shutil
import sys
import requests
import subprocess

from pathovar import get_external_databases_config, INSTALL_DIR

class SetupManager(object):
    def __init__(self, database_data, **opts):
        self.database_data = database_data
        self.verbose = opts.get('verbose', False)
        self.storage_path = self.database_data['storage_path']
        self.opts = opts

    def _pre_setup_dirs(self, *args, **kwargs):
        pass
    def _post_setup_dirs(self, *args, **kwargs):
        pass

    def setup_dirs(self, *args, **kwargs):
        if self.verbose: print("Setting up directories")
        try:
            os.makedirs(os.path.join(INSTALL_DIR, self.storage_path))
        except OSError, e:
            if self.verbose: print("Ignoring: ", e)
        self._pre_setup_dirs(*args, **kwargs)
        for key, value in self.database_data['data_urls'].items():
            if type(value) == list:
                try:
                    if self.verbose: print("Setting up %s" % os.path.join(INSTALL_DIR, self.storage_path, key))
                    os.makedirs(os.path.join(INSTALL_DIR, self.storage_path, key))
                except OSError, e:
                    if self.verbose: print("Ignoring: ", e)
        self._post_setup_dirs(*args, **kwargs)

    def _pre_download_file(self, *args, **kwargs):
        pass

    def _post_download_file(self, *args, **kwargs):
        pass

    def download_files(self, *args, **kwargs):
        for key, value in self.database_data['data_urls'].items():
            if type(value) == list:
                destination = os.path.join(INSTALL_DIR, self.storage_path, key)
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
        file_data.raise_for_status()
        file_name = url.split('/')[-1]
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
        self.setup_dirs(*args, **kwargs)
        self.download_files(*args, **kwargs)

    def remove(self, *args, **kwargs):
        shutil.rmtree(os.path.join(INSTALL_DIR, self.storage_path))
