import re
import json
from collections import defaultdict
from copy import copy

from bs4 import BeautifulSoup

HEADINGS = [
    "Pathways",
    "GeneOntology",
    "Interactions",
]

from pathovar.web.ncbi_xml import to_text_strip, to_int, to_attr_value

class GenBankGeneFile(object):
    CACHE_SCHEMA_VERSION = '0.3.5'
    def __init__(self, data, **opts):
        self.parser = None
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.gid = None
        self.gene_db_id = None
        self.gene_ref_tag = None
        self.gene_symbol = None
        self.org_name = None
        self.comments = None
        self.pathways = []
        self.go_terms = defaultdict(list)

        if "xml" in opts:
            self._parse_xml(data)
        elif "json" in opts:
            self._from_json(data)

    def _parse_xml(self, data):
        self.parser = BeautifulSoup(data)
        self.gene_db_id = self.parser.find("gene-track_geneid").get_text()
        self.gene_symbol = self.parser.find("gene-ref_locus")
        if self.gene_symbol:
            self.gene_symbol = to_text_strip(self.gene_symbol)
        self.gene_ref_tag = (self.parser.find("gene-ref_locus-tag"))
        if self.gene_ref_tag:
            self.gene_ref_tag = to_text_strip(self.gene_ref_tag)
        
        self.comments = [GenBankGeneFileComment(x.parent,self, xml=True) for x in self.parser.find_all('gene-commentary_heading') 
            if to_text_strip(x) in HEADINGS]
        for comment in self.comments:
            self.pathways.extend(comment.pathways)
            for key, value in comment.go_terms.items():
                self.go_terms[key].extend(value)

    def _from_json(self, json_dict):
        self.gene_db_id = json_dict['gene_db_id']
        self.gene_symbol = json_dict["gene_symbol"]
        self.gene_ref_tag =  json_dict["gene_ref_tag"]
        self.pathways = [Pathway(**data) for data in json_dict["pathways"]]
        self.go_terms = {category:[GOTerm(**data) for data in terms] for category, terms in  json_dict["go_terms"].items()}


    def to_json_safe_dict(self):
        data_dict = copy(self.__dict__)
        data_dict.pop("parser")
        data_dict.pop("comments")
        data_dict.pop("opts")
        data_dict['pathways'] = [pathway.to_json_safe_dict() for pathway in self.pathways]
        data_dict['go_terms'] = {category: [term.to_json_safe_dict() for term in terms] for category, terms in self.go_terms.items()}
        data_dict['schema_version'] = GenBankGeneFile.CACHE_SCHEMA_VERSION
        return data_dict

class GenBankGeneFileComment(object):
    def __init__(self, data, owner, **opts):
        self.owner = owner
        self.opts = opts
        self.verbose = opts.get('verbose', False)

        self.pathways = []
        self.go_terms = {}

        self._parse_xml(data)

    def _parse_xml(self, data):
        self.parser = data
        self.headings = [to_text_strip(h) for h in self.parser.find_all("gene-commentary_heading")]

        if "Pathways" in self.headings:
            self._parse_pathways()

        if "GeneOntology" in self.headings:
            self._parse_gene_ontologies()

    def _parse_pathways(self):
        pathways = []
        subcomponents = [x for x in self.parser.find_all("gene-commentary")]
        for subcomponent in subcomponents:
            text = to_text_strip(subcomponent.find("gene-commentary_text"))
            db = to_text_strip(subcomponent.find("dbtag_db"))
            id = to_text_strip(subcomponent.find("object-id_str"))
            url = to_text_strip(subcomponent.find("other-source_url"))
            pathway = Pathway(text, db, id, url)
            pathways.append(pathway)
        self.pathways = pathways

    def _parse_gene_ontologies(self):
        subcomponents = [x for x in self.parser.find_all("gene-commentary")]
        subcomponents_categories = [x.find("gene-commentary_label") for x in subcomponents if x.find("gene-commentary_label")]
        subcomponents_categories = {to_text_strip(x): x.parent for x in subcomponents_categories}

        go_terms = dict()
        for category_name in subcomponents_categories:
            category = subcomponents_categories[category_name]
            entries = category.find_all("gene-commentary")
            entry_terms = []
            for ent in entries:
                db = ent.find("dbtag_db").get_text()
                id = ent.find("object-id_id").get_text()
                term = ent.find("other-source_anchor").get_text()
                go_term = GOTerm(db, id, term)
                entry_terms.append(go_term)
            go_terms[category_name] = entry_terms
        
        self.go_terms = go_terms

    def __repr__(self):
        rep = "Comment(Headings: %(headings)s Content: %(pathways)r|%(go_terms)r)" % self.__dict__
        return rep

class Pathway(object):
    def __init__(self, name, db, id, url):
        self.name = name
        self.db = db
        self.id = id
        self.url = url
    
    def __repr__(self):
        rep = "Pathway(%(db)s|%(id)s|%(name)s)" % self.__dict__
        return rep

    def to_json_safe_dict(self):
        return self.__dict__

class GOTerm(object):
    def __init__(self, db, id, term):
        self.db = db
        self.id = id
        self.term = term

    def to_json_safe_dict(self):
        return self.__dict__

    def __repr__(self):
        rep = "GOTerm(%(db)s|%(id)s|%(term)s)" % self.__dict__
        return rep

class GenBankBioSystemFile(object):
    CACHE_SCHEMA_VERSION = "0.1b"
    def __init__(self, data, **opts):
        self.parser = None
        self.opts = opts
        self.system_names = None
        self.system_id = None
        self.system_description = None
        self.external_url = None
        self.external_accession = None
        self.system_categories = None

        if "xml" in self.opts:
            self._parse_xml(data)
        elif "json" in self.opts:
            self._from_json(data)

    def _parse_xml(self, data):
        self.parser = BeautifulSoup(data)
        self.system_names = [to_text_strip(name) for name in self.parser.find_all("system_names_e")]
        self.system_id = self.parser.find("sys-id_bsid").get_text()
        self.system_description = self.parser.find("system_description").get_text()
        self.external_url = self.parser.find("system_recordurl").get_text()
        self.external_accession = self.parser.find("system_externalaccn").get_text()
        self.system_categories = [to_text_strip(cat) for cat  in self.parser.find_all("system_category_e")]

    def to_json_safe_dict(self):
        data_dict = self.__dict__
        data_dict.pop("parser", None)
        
        data_dict['schema_version'] = GenBankBioSystemFile.CACHE_SCHEMA_VERSION
        return data_dict

    def _from_json(self, json_dict):
        self.system_names = json_dict["system_names"]
        self.system_id = json_dict["system_id"]
        self.system_description = json_dict["system_description"]
        self.external_url = json_dict["external_url"]
        self.external_accession = json_dict["external_accession"]
        self.system_categories = json_dict["system_categories"]

