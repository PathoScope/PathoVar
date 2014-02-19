##@package entrez_eutils
# Defines logic for interacting with Entrez's EUtils web services. Used to fetch genome 
# and sequence annotations. 
#
# TODO:
# - Migrate All components of the GenBankFeatureFile
# system to crossreference_snp_location.py in snp_annotation/
# - Include sequence data in GenBankFeature for first approximation
# of mutation prediction

# System Dependencies
import sys
import re
import os
from time import time, sleep
from collections import namedtuple

# External Dependencies
import requests
from bs4 import BeautifulSoup

from pathovar.web import get_robust

## URL CONSTANTS
taxonomy_detail_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={tid}'
# Organism name to Taxonomy ID: 
org_name_to_taxonomy_id = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={org_name}[organism]"

# Genome Link by Org Name
genome_by_org_name = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term={org_name}[organism]'

# Set remote environment for translating from genome to nucleotides
link_from_genome_to_nuccore =   'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id={id}&cmd=neighbor_history'
link_from_protein_to_gene =     "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=protein&db=gene&id={id}&cmd=neighbor_history"
link_from_protein_to_biosystem = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?&db=biosystems&dbfrom=protein&id={id}"
# Retrieve all matching nucleotide sequences
get_nucleotides_from_link = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&' \
                            'query_key={query_key}&WebEnv={web_env}&rettype={form}&retmode={mode}'

get_gene_from_link = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&rettype={form}&query_key={query_key}&WebEnv={web_env}"

# Retrieve a particular sequence record by gene id
get_nucleotides_by_gene_id = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={gene_id}&rettype={form}&retmode={mode}'

# Retrieve a particular sequence record by gene id
get_protein_by_gene_id = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={gene_id}&rettype={form}&retmode={mode}'

get_biosystem_by_bsid = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosystems&id={bsid}"
## EntrezEUtilsDriver
# Defines all of the logic for interacting with Entrez's EUtils web service, handling errors, 
# and performing multi-request actions.
#
class EntrezEUtilsDriver(object):
    def __init__(self, **opts):
        self.opts = opts
        self.verbose = self.opts.get('verbose', False)

    def find_genome_by_org_name(self, org_name, form = 'fasta', mode = 'text'):
        # Find the genome id
        if self.verbose: print("Fetching Genome ID from Entrez")
        genome_db_response = get_robust(genome_by_org_name.format(**dict(org_name = org_name)))
        genome_db_response.raise_for_status()
        
        # Parse the response
        if self.verbose: print("Response recieved...")
        genome_db_response_xml = BeautifulSoup(genome_db_response.text)
        genome_id = genome_db_response_xml.find_all('id')
        if(len(genome_id) == 0):
            if self.verbose: print("No Genome ID found")
            output_message = genome_db_response_xml.find_all('outputmessage')[0]
            if output_message.get_text() == u"No items found.":
                raise OrganismNotFoundException()
        # If there is ambiguity over which genome id to choose, just take the first
        genome_id = genome_id[0].get_text()
        if self.verbose: print("Genome ID: %s" % genome_id)

        # Set up remote environment to cross-link from genome id to nuccore id
        if self.verbose: print("Fetching Link to Genome from Entrez")
        genome_to_nuccore_response = get_robust(link_from_genome_to_nuccore.format(**dict(id = genome_id)))
        genome_to_nuccore_response.raise_for_status()

        if self.verbose: print("Response recieved...")
        query_key = None
        web_env = None        
        genome_to_nuccore_response_xml = BeautifulSoup(genome_to_nuccore_response.text)
        try:
            query_key = genome_to_nuccore_response_xml.find('querykey').get_text()
            web_env = genome_to_nuccore_response_xml.find('webenv').get_text()
        except:
            if(self.verbose): print(genome_to_nuccore_response_xml)
            raise EntrezEUtilsDriverException("Query Key and/or Web Env Missing")
        if self.verbose: print("Fetching Genome fasta file from Entrez")
        get_nucleotide_sequences_response = get_robust(get_nucleotides_from_link \
            .format(**dict(query_key = query_key, web_env = web_env, form = form, mode = mode)))
        get_nucleotide_sequences_response.raise_for_status()
        # The nucleotide sequence should be located in get_genome_response.text
        if self.verbose: print("Response recieved...")
        return get_nucleotide_sequences_response.text

    def find_nucleotides_by_gene_id(self, gene_id, form = 'fasta', mode = 'text'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        timer = time()
        if self.verbose: print("Fetching %s from Entrez" % gene_id)
        get_nucleotide_sequences_response = get_robust(get_nucleotides_by_gene_id.format(**dict(gene_id=gene_id, form=form, mode=mode)))
        if self.verbose: print("Response recieved (%s sec)" % str(time() - timer))
        get_nucleotide_sequences_response.raise_for_status()
        # The genome sequence should be located in get_genome_response.text
        return get_nucleotide_sequences_response.text

    def find_protein_by_gene_id(self, gene_id, form = 'fasta', mode = 'text'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        timer = time()
        if self.verbose: print("Fetching %s from Entrez" % gene_id)
        get_protein_sequences_response = get_robust(get_protein_by_gene_id.format(**dict(gene_id=gene_id, form=form, mode = mode)))
        if self.verbose: print("Response recieved (%s sec)" % str(time() - timer))
        get_protein_sequences_response.raise_for_status()
        return get_protein_sequences_response.text
        
    def find_gene_by_gene_id(self, gene_id, form = 'xml'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        timer = time()
        if self.verbose: print("Fetching %s from Entrez" % gene_id)
        get_link_from_protein_to_gene_response = get_robust(link_from_protein_to_gene.format(**dict(id=gene_id)))
        if self.verbose: print("Response recieved (%s sec)" % str(time() - timer))
        get_link_from_protein_to_gene_response.raise_for_status()
        query_key = None
        web_env = None
        get_link_from_protein_to_gene_response_xml = BeautifulSoup(get_link_from_protein_to_gene_response.text)
        try:
            query_key = get_link_from_protein_to_gene_response_xml.find('querykey').get_text()
            web_env = get_link_from_protein_to_gene_response_xml.find('webenv').get_text()
        except:
            if(self.verbose): print(get_link_from_protein_to_gene_response_xml)
            raise EntrezEUtilsDriverException("Query Key and/or Web Env Missing")
        get_gene_from_link_response = get_robust(get_gene_from_link.format(**dict(id=gene_id, form=form, query_key=query_key, 
            web_env=web_env)))
        get_gene_from_link_response.raise_for_status()
        return get_gene_from_link_response.text

    def find_biosystem_ids_by_gene_id(self, gene_id):
        gene_id = str(gene_id)
        timer = time()
        #if self.verbose: print("Fetching %s's BioSystems from Entrez" % gene_id)
        link_from_protein_to_biosystem_response = get_robust(link_from_protein_to_biosystem.format(**dict(id=gene_id)))
        link_from_protein_to_biosystem_response.raise_for_status()
        link_set = BeautifulSoup(link_from_protein_to_biosystem_response.text)
        link_set_ids = set([link.get_text() for link in link_set.find_all("id") if link.get_text() != gene_id])
        return link_set_ids

    def find_biosystem_by_bsid(self, bsid):
        bsid = str(bsid)
        timer = time()
        if self.verbose: print("Fetching %s from Entrez" % bsid)
        get_biosystem_by_bsid_response = get_robust(get_biosystem_by_bsid.format(**dict(bsid=bsid)))
        get_biosystem_by_bsid_response.raise_for_status()
        return get_biosystem_by_bsid_response.text

## EntrezEUtilsDriverException
# Parent class for capturing all EUtils generated exceptions. Exception class 
# representing programmatic errors while fetching information with EntrezEUtilsDriver
class EntrezEUtilsDriverException(Exception):
    pass

## OrganismNotFoundException
# Exception indicating the organism indicated by org_name was not 
# found on Entrez
class OrganismNotFoundException(EntrezEUtilsDriverException):
    pass