import glob
import os
import re
import subprocess
import sys
import tempfile

from collections import defaultdict

from pathovar.web import entrez_eutils
from pathovar.utils import defline_parser

def main(snpPath, vcfFile, tempDir = None, gidMap = None, **opts):
	#JOSH: modify the following 4 lines to instead be passed as arguments to the method
	#Can rename the method to whatever you want if main won't work
	#Please use gidMap = {} as the default if the argument isn't given
	if tempDir == None:
		tempDir = tempfile.mkdtemp()
	if gidMap == None:
		gidMap = dict()
	#Check  snpEff permissions
	snpFile = snpPath + "/snpEff.jar"
	snpData = snpPath + "/data/"
	verbose = opts.get("verbose", False)
	if not os.path.exists(tempDir):
		raise snpEffConfigException("ERROR:Output directory invalid/no permissions")
	if not os.path.exists(snpPath):
		raise snpEffConfigException("ERROR:Invalid snpEff path/permissions")
	if not os.path.isfile(snpFile):
		raise snpEffConfigException("ERROR:snpEff.jar does not exist in specified directory")
	streamMap= {}
	#separate original into vcfs into sub vcfs
	try:
		with open(vcfFile,"r") as inStream:
			header = ""
			buildHead = 1
			dumpStream = open(tempDir + "/failedlines_pathovar.vcf",'w')
			for line in inStream:
				if buildHead and re.search(r'^#',line):
					header = header + line
				elif re.search(r'^[ti,gi][^\t]*\t',line):
					buildHead = 0
					parseObj = defline_parser(line.split("\t")[0])
					if (not parseObj.has_key("gene_id")):
						dumpStream.write(line)
						continue
					gID = parseObj["gene_id"]
					strainName = ""
					#manually look for strain name by index if needed
					if not gidMap.has_key(gID):
						tag = line.split()[0]
						splitTag = tag.split("|")
						if len(splitTag) ==  9:
							gidMap[gID] = splitTag[7]
					strainName = gidMap[gID]
					if gidMap.has_key(gID):
						modLine = re.sub(r'^[ti,gi][^\t]*\t',strainName.split(".")[0] + "\t",line)
						#create stream and file if needed
						if not (gID in streamMap.keys()):
							addStream = open(tempDir + "/" +str(gID)+"_pathovar.vcf",'w')
							#addStream.write(header)
							streamMap[gID] = addStream
						writeStream = streamMap[gID]
						writeStream.write(modLine)
					else:
						dumpStream.write(line)
			dumpStream.close()
			for streamOpen in streamMap.values():
				streamOpen.close()
	except Exception, e:
		print "Exception occurred during file separation"
		print e
	#build databases for all files
	buildErrorSet = set()
	validSet = set()
	for gID in streamMap.keys():
		try:
			buildDatabase(snpPath,gID,gidMap[gID])	
			validSet.add(gID)
		except Exception, e:
			print e
			buildErrorSet.add(gID)
	annoErrorSet = set()
	finalSet = {}
	for geneID in validSet:
		try:
			annoFile = annotate(snpPath,gidMap[geneID],streamMap[geneID].name)
			finalSet[geneID] = annoFile
		except Exception,e:
			print e
			annoErrorSet.add(geneID)
	print(buildErrorSet)
	print(annoErrorSet)
	resultsDict = defaultdict(lambda : defaultdict(int))
	for geneID in finalSet.keys():
		try:
			with open(finalSet[geneID],"r") as curFile:
				for line in curFile:
					if not re.search(r'^#',line):
						splitArr = line.split()
						loc = splitArr[1]
						effArr = splitArr[7].split(";EFF=")
						if len(effArr) > 1:
							resultsDict[geneID][loc] = parseTag(effArr[1])
		except Exception:
			pass

	cleanDir(tempDir)
	return resultsDict

#Call to get the necessary database files
#THROWS exceptions from entrez eutils
#snpPath is path to snpEff directory
#geneID is the gid of the strain
#strainName is the accession name
def fetchGenome(snpPath,genomeDir,geneID,strainName):
	eutils_handle = entrez_eutils.EntrezEUtilsDriver()
	if not os.path.isfile(genomeDir + "/" + strainName + '.fa'):
		genomeSequence = eutils_handle.find_nucleotides_by_gene_id(geneID)
		splitSeq = genomeSequence.split("\n")
		splitSeq[0] = ">" + strainName
		writeStream = open(genomeDir +"/"+ strainName + '.fa', 'w')
		for line in splitSeq:
			writeStream.write(line+"\n")
	if not os.path.isfile(genomeDir + '/genes.gb'):
		genomeAnnotation = eutils_handle.find_nucleotides_by_gene_id(geneID,mode='gb',form='genbank')
		open(genomeDir +"/"+ 'genes.gb','w').write(genomeAnnotation)
	needsMod = 1
	with open(snpPath+"/snpEff.config", "r") as configfile:
		for line in configfile:
			if(strainName + ".genome")  in line:
				needsMod = 0
	if needsMod:
		addFile = open(snpPath+"/snpEff.config","a")
		addFile.write("%s.genome : %s\n" % (strainName,strainName))
		addFile.close()

#Call to build a database. Will download file if needed
#THROWS subprocess.CalledProcessError if building doesn't work
#snpPath is path to snpEff directory
#geneID is the gid of the strain
#strainName is the accession name
def buildDatabase(snpPath,geneID,strainName):
	genomeDir = snpPath + "/data/" + strainName + "/"
	snpFile = snpPath + "/snpEff.jar"
	if not os.path.exists(genomeDir):
		try:
			os.makedirs(genomeDir)
		except OSError:
			raise snpEffPermissionsException("WARNING: Unable to create directory in snpEff due to permissions")
	if not os.path.isfile(genomeDir + '/snpEffectPredictor.bin'):
		fetchGenome(snpPath,genomeDir,geneID,strainName)
		subprocess.check_call(["java","-jar",snpFile,"build","-genbank",strainName])

#Call to annotate a file. SHOULD BUILD DATABASE FIRST
#snpPath is path to snpEff directory
#strainName is the accession name
#splitFile is the .vcf file containing only the strains SNP's and chromosomes
#removeFlag signals to delete the .vcf file after annotating it
def annotate(snpPath,strainName,splitFile):
	snpFile = snpPath + "/snpEff.jar"
	outFile = splitFile.split(".vcf")[0] + ".eff.vcf"
	cmdLine = "java -Xmx2g -jar " + snpFile + " -c " + snpPath + "/snpEff.config " + strainName + " -noStats -no-downstream -no-intron -no-upstream " + splitFile + " > " + outFile
	subprocess.check_call(cmdLine,shell=True)
	return outFile

#Call to parse a snpEff tag into components
#returns a dict of dict object
#dictObj[changeNumber][info] is the format of the returned dict
#dictObj.keys() will extract the keys for each of the changes
#dictObj[changeNum]["type"] is the type of mutation
#dictObj[changeNum]["gene"] is the gene the mutation was in, returns 0 if no gene
#dictObj[changeNum]["error"] is the type of error, returns 0 if no errors
def parseTag(extractTag):
	annotations = extractTag.split(",")
	infoReturn = defaultdict(lambda: defaultdict(int))
	nextKey = 0
	for change in annotations:
		infoReturn[nextKey]["type"] = change.split("(")[0]
		compArr = change.split("|")
		gene = compArr[5]
		aaChange = compArr[3]
		if len(gene) > 0:
			infoReturn[nextKey]["gene"] = gene
		if len(compArr) > 11:
			infoReturn[nextKey]["error"] = compArr[11].split(")")[0]
		if len(aaChange) > 0:
			infoReturn[nextKey]["aaChange"] = aaChange
		nextKey = nextKey + 1
	return infoReturn

def cleanDir(dir):
	map(os.remove, glob.glob(os.path.join(dir,"*.vcf")))
	os.rmdir(dir)

#Cant find the config file
class snpEffConfigException(Exception):
	pass

#Cant run the methods
class snpEffRunException(Exception):
	pass

#Cant write to directory
class snpEffPermissionsException(Exception):
	pass

if __name__ == "__main__":
	main()
