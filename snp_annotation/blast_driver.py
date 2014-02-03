##@blast_driver
# Utilities for programmatically calling NCBI Blast+

import os
import sys
import subprocess

from bs4 import BeautifulSoup

from pathovar.utils import defline_parser

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

NUCLEOTIDE = "nucl"
PROTEIN = "prot"
class BlastAnnotationDriver(object):

	def __init__(self, database_path, database_type, bin_dir = '', **opts):
		self.database_path = database_path
		self.database_type = database_type
		self.bin_dir = bin_dir
		self.opts = opts
		self.verbose = opts.get('verbose', False)
		self.process = None
		self.db_name = os.path.splitext(os.path.basename(self.database_path))[0]
		self.outfile = None
		if not self.is_built():
			self.build_database()

	def results(self):
		if self.process.returncode == 0 and outfile:
			return BlastResultsXMLParser(self.outfile)
		else:
			self.process.communicate()
			return self.results()

	def is_built(self):
		if self.database_type == NUCLEOTIDE:
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
		if not self.is_built():
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


class NucleotideDatabaseBlastAnnotatorBase(object):
	def __init__(self, database_file_paths, **opts):
		self.opts = opts
		self.verbose = opts.get('verbose', False)
		self.blast_drivers = map(lambda dbf: BlastAnnotationDriver(dbf, NUCLEOTIDE), database_file_paths)

	@property
	def processes(self):
		return [b.process for b in self.blast_drivers]

	def wait_for_results(self):
		ret_codes = [b.process.wait() for b in self.blast_drivers]
		if sum(ret_codes) != 0:
			raise BlastDriverException("A %s Blast job failed %r" % (str(type(self)), ret_codes))
		result_files_dict = {self.blast_drivers[ind].db_name:BlastResultsXMLParser(outfile) for ind, outfile in enumerate(self.outfiles)}
		return result_files_dict

	@property
	def outfiles(self):
		return [b.outfile for b in self.blast_drivers]

	def query_with_nucleotides(self, query, **opts):
		for blaster in self.blast_drivers:
			if self.verbose: print("Blasting against %s" % blaster.db_name)
			proc = blaster.blast_with_nucleotides(query, **opts)

class ProteinDatabaseBlastAnnotatorBase(object):
	def __init__(self, database_file_paths, **opts):
		self.opts = opts
		self.verbose = opts.get('verbose', False)
		self.blast_drivers = map(lambda dbf: BlastAnnotationDriver(dbf, PROTEIN), database_file_paths)

	@property
	def processes(self):
		return [b.process for b in self.blast_drivers]

	def wait_for_results(self):
		ret_codes = [b.process.wait() for b in self.blast_drivers]
		if sum(ret_codes) != 0:
			raise BlastDriverException("A %s Blast job failed %r" % (str(type(self)), ret_codes))
		result_files_dict = {self.blast_drivers[ind].db_name:BlastResultsXMLParser(outfile) for ind, outfile in enumerate(self.outfiles)}
		return result_files_dict

	@property
	def outfiles(self):
		return [b.outfile for b in self.blast_drivers]

	def query_with_proteins(self, query, **opts):
		for blaster in self.blast_drivers:
			if self.verbose: print("Blasting against %s" % blaster.db_name)
			proc = blaster.blast_with_proteins(query, **opts)

## Blast Results Parsing
class BlastResultsXMLParser(object):
	def __init__(self, file_path, **opts):
		self.file_path = file_path
		self.opts = opts
		self.parser = BeautifulSoup(''.join(open(file_path).readlines()))
		self.queries = [BlastResultsQuery(iteration) for iteration in self.parser.find_all('iteration')]
		self.queries = {defline_parser(query.query_def)['gene_id']: query for query in self.queries if len(query.hits) > 0}


class BlastResultsQuery(object):
	def __init__(self, parser):
		self.parser = parser
		self.query_def = parser.find("iteration_query-def").get_text()
		self.__len__ = len(parser.find("iteration_query-len").get_text())
		self.hits = sorted([BlastResultsHit(h) for h in parser.find_all("hit")], key=lambda x: x.e_value)

	def __repr__(self):
		rep = 'BlastResultQuery(%s Hits:%d)' % (self.query_def, len(self.hits))
		return rep

class BlastResultsHit(object):
	def __init__(self, parser):
		self.parser = parser
		self.hit_def = parser.find('hit_def').get_text()
		self.__len__ = int(parser.find("hit_len").get_text())
		self.hsps = [BlastResultsHSP(hsp) for hsp in parser.find_all("hsp")]
		self.e_value = min(self.hsps, key = lambda x: x.e_value).e_value

	def __repr__(self):
		rep = "BlastResultHit(%s E Value:%0.000f HSPs:%d)" % (self.hit_def, self.e_value, len(self.hsps))
		return rep

class BlastResultsHSP(object):
	def __init__(self, parser):
		self.parser = parser
		self.e_value = float(parser.find("hsp_evalue").get_text())
		# Query range
		self.query_from = int(parser.find("hsp_query-from").get_text())
		self.query_to = int(parser.find("hsp_query-to").get_text())

		# Hit range
		self.hit_from = int(parser.find("hsp_hit-from").get_text())
		self.hit_to = int(parser.find("hsp_hit-to").get_text())

		# Frames (Nucleotide Derived Only)
		try:
			self.query_frame = parser.find('hsp_query-frame').get_text()
		except ValueError, e:
			pass
		try:
			self.hit_frame = parser.find('hsp_hit-frame').get_text()
		except ValueError, e:
			pass

		# Alignment
		self.identity = parser.find("hsp_identity").get_text()
		self.gaps = int(parser.find('hsp_gaps').get_text())
		self.__len__ = int(parser.find('hsp_align-len').get_text())
		self.qseq = parser.find("hsp_qseq").get_text()
		self.hseq = parser.find("hsp_hseq").get_text()
		self.midline = parser.find("hsp_midline").get_text()

	def __repr__(self):
		rep = "BlastResultsHSP(E Value:%0.000f)" % self.e_value
		return rep

class BlastDriverException(Exception):
	pass



if __name__ == '__main__':
	
	main(sys.argv)