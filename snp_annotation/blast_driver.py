##@blast_driver
# Utilities for programmatically calling NCBI Blast+

import os
import sys
import subprocess

## Location of the ncbi-blast+ binaries
BLAST_BIN_DIR = ""

def makeblastdb(in_file = None, dbtype='nucl'):
	if in_file == None:
		raise IOError("makeblastdb -in file not found")
	args = {"in": in_file, "dbtype": dbtype}
	subprocess.call(BLAST_BIN_DIR + "makeblastdb" + ' -in {in} -dbtype {dbtype}'.format(**args), shell = True)

def blastn(query, db_name, evalue = 0.001, num_threads = 1, outfmt = 5, outfile = None):
	if query == None or not os.path.exists(query):
		raise IOError("blastn -query file not found")
	if db_name == None or not os.path.exists(db_name):
		raise IOError("blastn -db_name file not found")
	if outfile == None:  
		outfile = query + "_vs_" + db_name + '.blastn'
	args = {"query" : query, "db" : db_name, "evalue" : evalue, "num_threads": num_threads, "outfmt": outfmt, 'out': outfile}

	call = subprocess.Popen(BLAST_BIN_DIR + "blastn -query {query} -db {db} -evalue {evalue} -num_threads {num_threads} -outfmt {outfmt} -out {out}".format(**args), 
		stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
	return call

def main(args):
	#makeblastdb(args[1])
	stdout, stderr = blastn(args[1], args[2]).communicate()
	print(stdout)
	print(stderr)


if __name__ == '__main__':
	main(sys.argv)