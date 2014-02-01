import os
import sys
import requests
import subprocess

from pathovar import get_external_databases_config, INSTALL_DIR

database_data = get_external_databases_config()['drugbank']

def setup(data_urls, destination_dir):
    try:
        os.makedirs(os.path.join(INSTALL_DIR, destination_dir))
    except:
        pass

    for url in data_urls:
        setup_file(url, destination_dir)

def setup_file(url, destination_dir):
    print("Getting %s" % url)
    file_data = requests.get(url)
    file_data.raise_for_status()
    file_name = url.split('/')[-1]
    with open(destination_dir + os.sep + file_name, 'wb') as data_file:
        data_file.write(file_data.content)
    unzip = os.path.splitext(file_name)
    if unzip[1] == '.zip':
        os.system('unzip ' + destination_dir + os.sep + file_name + ' -d' + destination_dir)
    return destination_dir + os.sep + file_name[:-4]

if __name__ == '__main__':
    dest = sys.argv[1]
    data_urls = sys.argv[2:]
    setup(data_urls, dest)
