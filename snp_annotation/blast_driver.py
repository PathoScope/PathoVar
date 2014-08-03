##@blast_driver
# Utilities for programmatically calling NCBI Blast+

import os
import sys
import subprocess
from copy import copy

from pathovar.web.ncbi_xml import ET

from pathovar.utils import defline_parser

NUCLEOTIDE = "nucl"
PROTEIN = "prot"
class BlastAnnotationDriver(object):

	def __init__(self, database_path, database_type, bin_dir = '', **opts):
		self.database_path = database_path
		self.database_type = database_type
		self.bin_dir = bin_dir
		if self.bin_dir != "":
			self.bin_dir += os.sep
		self.opts = opts
		self.verbose = opts.get('verbose', False)
		self.process = None
		self.db_name = os.path.splitext(os.path.basename(self.database_path))[0]
		self.outfile = None
		if not self.is_built():
			self.build_database()

	def clean(self):
		if self.opts.get("clean", False):
			if self.process.returncode == 0 and self.outfile:
				print("Removing %s" % self.outfile)
				os.remove(self.outfile)
				self.outfile = None

	def results(self):
		if self.process.returncode == 0 and self.outfile:
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
		if not self.is_built() or proc != 0:
			raise IOError("Blast Database: makeblastdb Failure (%s). Is NCBI-Blast+ on your path or in your configuration file?" % self.database_path)

	def _run_blast(self, query_file, cmd, **opts):
		if not self.is_built():
			raise BlastDriverException("Database Not Built. Cannot run %s on %s" % (cmd, self.database_path))
		args = dict(evalue = "0." + ("0" * 50) + "1", num_threads = 2, outfmt = 5, 
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
	def __init__(self, database_file_paths, db_name_prefix = '', bin_dir = '',**opts):
		self.opts = opts
		self.verbose = opts.get('verbose', False)
		if(len(database_file_paths) == 0):
			raise BlastDriverException("No database files were found!")
		self.blast_drivers = map(lambda dbf: BlastAnnotationDriver(dbf, NUCLEOTIDE, bin_dir, **opts), database_file_paths)
		self.collection_name = db_name_prefix
		for driver in self.blast_drivers:
			driver.db_name = db_name_prefix + '_' + driver.db_name

	@property
	def processes(self):
		return [b.process for b in self.blast_drivers]

	def wait_for_results(self):
		ret_codes = [b.process.wait() for b in self.blast_drivers]
		if sum(ret_codes) != 0:
			if self.verbose:
				print("BlastDriver Error")
				for i, rc in enumerate(ret_codes):
					if rc != 0:
						errd_blast_driver = self.blast_drivers[i]
						print("Driver %d" % i)
						print("STDOUT")
						print(errd_blast_driver.process.stdout.readlines())
						print("STDERR")
						print(errd_blast_driver.process.stderr.readlines())
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

	def clean(self):
		for blaster in self.blast_drivers:
			blaster.clean()


class ProteinDatabaseBlastAnnotatorBase(object):
	def __init__(self, database_file_paths, db_name_prefix = '', bin_dir = '',**opts):
		self.opts = opts
		self.verbose = opts.get('verbose', False)
		if(len(database_file_paths) == 0):
			raise BlastDriverException("No database files were found!")
		self.blast_drivers = map(lambda dbf: BlastAnnotationDriver(dbf, PROTEIN, bin_dir, **opts), database_file_paths)
		self.collection_name = db_name_prefix
		for driver in self.blast_drivers:
			driver.db_name = db_name_prefix + '_' + driver.db_name

	@property
	def processes(self):
		return [b.process for b in self.blast_drivers]

	def wait_for_results(self):
		ret_codes = [b.process.wait() for b in self.blast_drivers]
		if sum(ret_codes) != 0:
			if self.verbose:
				print("BlastDriver Error")
				for i, rc in enumerate(ret_codes):
					if rc != 0:
						errd_blast_driver = self.blast_drivers[i]
						print("Driver %d" % i)
						print("STDOUT")
						print(errd_blast_driver.process.stdout.readlines())
						print("STDERR")
						print(errd_blast_driver.process.stderr.readlines())

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

	def clean(self):
		for blaster in self.blast_drivers:
			blaster.clean()

## Blast Results Parsing
# Parses the XML output of the BLAST+ program for simple searches like blastn, blastp, or blastx.
# The resulting object is a 
class BlastResultsXMLParser(object):
	def __init__(self, file_path, **opts):
		self.file_path = file_path
		self.opts = opts
		self.parser = ET.fromstring(''.join(open(file_path).readlines()))
		self.queries = [BlastResultsQuery(iteration) for iteration in self.parser.findall('.//Iteration')]
		self.queries = {defline_parser(query.query_def)['gene_id']: query for query in self.queries if len(query.hits) > 0}

	def __iter__(self):
		for q, r in self.queries.iteritems():
			yield (q, r)

	def __len__(self):
		return len(self.queries)

	def to_json_safe_dict(self):
		return {k:v.to_json_safe_dict() for k, v in self.queries.items()}

class BlastResultsQuery(object):
	def __init__(self, parser):
		self.query_def = parser.find(".//Iteration_query-def").text
		self.__len__ = len(parser.find(".//Iteration_query-len").text)
		self.hits = sorted([BlastResultsHit(h) for h in parser.findall(".//Hit")], key=lambda x: x.e_value)

	def to_json_safe_dict(self):
		data_dict = copy(self.__dict__)
		data_dict['hits'] = [h.to_json_safe_dict() for h in self.hits]
		return data_dict

	def __repr__(self):
		rep = 'BlastResultQuery(%s Hits:%d)' % (self.query_def, len(self.hits))
		return rep

class BlastResultsHit(object):
	def __init__(self, parser):
		self.hit_def = parser.find('.//Hit_def').text
		self.__len__ = int(parser.find(".//Hit_len").text)
		self.hsps = [BlastResultsHSP(hsp) for hsp in parser.findall(".//Hsp")]
		self.e_value = min(self.hsps, key = lambda x: x.e_value).e_value

	def to_json_safe_dict(self):
		data_dict = copy(self.__dict__)
		data_dict['hsps'] = [hsp.to_json_safe_dict() for hsp in self.hsps]
		return data_dict


	def __repr__(self):
		rep = "BlastResultHit(%s E Value:%0.000f HSPs:%d)" % (self.hit_def, self.e_value, len(self.hsps))
		return rep

class BlastResultsHSP(object):
	def __init__(self, parser):
		
		self.e_value = float(parser.find(".//Hsp_evalue").text)
		# Query range
		self.query_from = int(parser.find(".//Hsp_query-from").text)
		self.query_to = int(parser.find(".//Hsp_query-to").text)

		# Hit range
		self.hit_from = int(parser.find(".//Hsp_hit-from").text)
		self.hit_to = int(parser.find(".//Hsp_hit-to").text)

		# Frames (Nucleotide Derived Only)
		try:
			self.query_frame = parser.find('.//Hsp_query-frame').text
		except ValueError, e:
			pass
		try:
			self.hit_frame = parser.find('.//Hsp_hit-frame').text
		except ValueError, e:
			pass

		# Alignment
		self.identity = parser.find(".//Hsp_identity").text
		self.gaps = int(parser.find('.//Hsp_gaps').text)
		self.__len__ = int(parser.find('.//Hsp_align-len').text)
		self.qseq = parser.find(".//Hsp_qseq").text
		self.hseq = parser.find(".//Hsp_hseq").text
		self.midline = parser.find(".//Hsp_midline").text

	def to_json_safe_dict(self):
		return self.__dict__

	def __repr__(self):
		rep = "BlastResultsHSP(E Value:%0.000f)" % self.e_value
		return rep

class BlastDriverException(Exception):
	pass
