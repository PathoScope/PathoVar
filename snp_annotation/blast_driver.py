##@blast_driver
# Utilities for programmatically calling NCBI Blast+

import os
import sys
import subprocess
from bs4 import BeautifulSoup

## Location of the ncbi-blast+ binaries
BLAST_BIN_DIR = ""

def make_blastdb(in_file = None, dbtype='nucl'):
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
	#make_blastdb(args[1])
	stdout, stderr = blastn(args[1], args[2]).communicate()
	print(stdout)
	print(stderr)

class BlastAnnotationDriver(object):
	NUCLEOTIDE = "nucleotide"
	PROTEIN = "protein"
	def __init__(self, database_path, database_type, bin_dir = '', **opts):
		self.database_path = database_path
		self.database_type = database_type
		self.bin_dir = bin_dir
		self.opts = opts
		self.verbose = opts.get('verbose', False)
		self.process = None
		self.db_name = os.path.splitext(os.path.basename(self.database_path))[0]
		self.outfile = None

	def results(self):
		if self.process.returncode == 0 and outfile:
			return BlastResultsXMLParser(self.outfile)
		else:
			self.process.communicate()
			return self.results()

	def is_built(self):
		if self.database_type == 'nucleotide':
			return (os.path.exists(self.database_path + '.nhr')	and
				os.path.exists(self.database_path + '.nsq') and
				os.path.exists(self.database_path + '.nin'))
		else:
			return (os.path.exists(self.database_path + '.phr')	and
				os.path.exists(self.database_path + '.psq') and
				os.path.exists(self.database_path + '.pin'))

	def build_database(self):
		if not os.path.exists(self.database_path):
			raise IOError("Blast Database (%s) could not be found!" % self.database_path)
		args = {'in': self.database_path, 'dbtype': self.database_type}
		proc = subprocess.call(self.bin_dir + "makeblastdb" + 
			' -in {in} -dbtype {dbtype}'.format(**args), shell = True)
		if proc != 0: 
			raise IOError("Blast Database: makeblastdb Failure (%s)" % self.database_path)

	def _run_blast(self, query_file, cmd, **opts):
		if not self.is_built():
			raise BlastDriverException("Database Not Built. Cannot run %s on %s" % (cmd, self.database_path))
		args = dict(evalue = 0.001, num_threads = 1, outfmt = 5, 
			outfile = query_file + "_on_" + self.db_name + "." + cmd + '.xml')
		for opt in opts:
			args[opt] = opts[opt]
		args['cmd'] = cmd
		args['query'] = query_file
		args['db'] = self.database_path
		self.outfile = args['outfile']
		self.process = subprocess.Popen(self.bin_dir + "{cmd} -query {query} -db {db} -evalue {evalue} -num_threads {num_threads} -outfmt {outfmt} -out {outfile}".format(**args), 
		stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
		return self.process

	def blast_with_nucleotides(self, query, **opts):
		cmd = 'blastn'
		if self.database_type == "protein":
			cmd = 'blastx'
		return self._run_blast(query, cmd, **opts)

	def blast_with_proteins(self, query, **opts):
		cmd = 'blastp'
		if self.database_type == "nucleotide":
			cmd = "tblastn"
		return self._run_blast(query, cmd, **opts)



class BlastResultsXMLParser(object):
	def __init__(self, file_path, **opts):
		self.file_path = file_path
		self.opts = opts
		self.parser = BeautifulSoup(''.join(open(file_path).readlines()))
		#self.hits = self.parser.find_all()


class BlastResultHit(object):
	pass

class BlastResulHSP(object):
	pass










class BlastDriverException(Exception):
	pass



if __name__ == '__main__':
	
	main(sys.argv)

