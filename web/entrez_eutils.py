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
from collections import namedtuple

# External Dependencies
import requests
from bs4 import BeautifulSoup

## URL CONSTANTS

taxonomy_summary_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={tid}'
taxonomy_detail_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={tid}'
# Organism name to Taxonomy ID: 
org_name_to_taxonomy_id = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term={org_name}[organism]"

# Genome Link by Org Name
genome_by_org_name = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=genome&term={org_name}[organism]'

# Set remote environment for translating from genome to nucleotides
link_from_genome_to_nuccore = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=genome&db=nuccore&id={id}&cmd=neighbor_history'

# Retrieve all matching nucleotide sequences
get_nucleotides_from_link = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&' \
                            'query_key={query_key}&WebEnv={web_env}&rettype={form}&retmode={mode}'
# Retrieve a particular sequence record by gene id
get_nucleotides_by_gene_id = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={gene_id}&rettype={form}&retmode={mode}'

# Retrieve a particular sequence record by gene id
get_protein_by_gene_id = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={gene_id}&rettype={form}&retmode={mode}'

get_gene_by_gene_id = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id={gene_id}&rettype={form}"

## EntrezEUtilsDriver
# Defines all of the logic for interacting with Entrez's EUtils web service, handling errors, 
# and performing multi-request actions.
#
class EntrezEUtilsDriver(object):
    def __init__(self, opts):
        self.opts = opts
        self.verbose = self.opts['verbose']

    def find_genome_by_org_name(self, org_name, form = 'fasta', mode = 'text'):
        # Find the genome id
        if self.verbose: print("Fetching Genome ID from Entrez")
        genome_db_response = requests.get(genome_by_org_name.format(**dict(org_name = org_name)))
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
        genome_to_nuccore_response = requests.get(link_from_genome_to_nuccore.format(**dict(id = genome_id)))
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
        get_nucleotide_sequences_response = requests.get(get_nucleotides_from_link \
            .format(**dict(query_key = query_key, web_env = web_env, form = form, mode = mode)))
        get_nucleotide_sequences_response.raise_for_status()
        # The nucleotide sequence should be located in get_genome_response.text
        if self.verbose: print("Response recieved...")
        return get_nucleotide_sequences_response.text

    def find_nucleotides_by_gene_id(self, gene_id, form = 'fasta', mode = 'text'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        if self.verbose: print("Fetching data from Entrez")
        get_nucleotide_sequences_response = requests.get(get_nucleotides_by_gene_id.format(**dict(gene_id=gene_id, form=form, mode=mode)))
        get_nucleotide_sequences_response.raise_for_status()
        # The genome sequence should be located in get_genome_response.text
        return get_nucleotide_sequences_response.text

    def find_protein_by_gene_id(self, gene_id, form = 'fasta', mode = 'text'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        if self.verbose: print("Fetching data from Entrez")
        get_protein_sequences_response = requests.get(get_protein_by_gene_id.format(**dict(gene_id=gene_id, form=form, mode = mode)))
        get_protein_sequences_response.raise_for_status()
        return get_protein_sequences_response.text
        
    def find_gene_by_gene_id(self, gene_id, form = 'xml'):
        # Just in case it was passed as an int
        gene_id = str(gene_id)
        if self.verbose: print("Fetching data from Entrez")
        get_gene_response = requests.get(get_gene_by_gene_id.format(**dict(gene_id=gene_id, form=form)))
        get_gene_response.raise_for_status()
        if self.verbose:
            open(gene_id+".gene." + form, 'w').write(get_gene_response.text)
        return get_gene_response.text

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