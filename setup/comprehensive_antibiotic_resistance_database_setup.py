import os
import sys
import requests
import subprocess

from pathovar.snp_annotation import blast_driver

def setup(data_urls, destination_dir):
    try:
        os.makedirs('databases/comprehensive_antibiotic_resistance_database')
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
    if unzip[1] == '.gz':
        os.system('gunzip ' + destination_dir + file_name)
        blast_driver.make_blastdb(destination_dir + file_name[:-3])
    return destination_dir + os.sep + file_name

if __name__ == '__main__':
    dest = sys.argv[1]
    data_urls = sys.argv[2:]
    setup(data_urls, dest)
