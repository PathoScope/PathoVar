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

CACHE_SCHEMA_VERSION = '0.3.4'

from pathovar.web.ncbi_xml import to_text_strip, to_int, to_attr_value

class GenbankGeneFile(object):
    def __init__(self, data, **opts):
        self.parser = None
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.gid = None
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
        self.gene_symbol = self.parser.find("gene-ref_locus")
        if self.gene_symbol:
            self.gene_symbol = self.gene_symbol.to_text_strip()
        self.gene_ref_tag = (self.parser.find("gene-ref_locus-tag"))
        if self.gene_ref_tag:
            self.gene_ref_tag = to_text_strip(self.gene_ref_tag)
        
        self.comments = [GenbankGeneFileComment(x.parent,self, xml=True) for x in self.parser.find_all('gene-commentary_heading') 
            if to_text_strip(x) in HEADINGS]
        for comment in self.comments:
            self.pathways.extend(comment.pathways)
            for key, value in comment.go_terms.items():
                self.go_terms[key].extend(value)

    def _from_json(self, json_dict):
        self.gene_symbol = json_dict["gene_symbol"]
        self.gene_ref_tag =  json_dict["gene_ref_tag"]
        self.pathways = [Pathway(**data) for data in json_dict["pathways"]]
        self.go_terms = {category:[GOTerm(**data) for data in terms] for category, terms in  json_dict["go_terms"].items}


    def to_json_safe_dict(self):
        data_dict = copy(self.__dict__)
        data_dict.pop("parser")
        data_dict.pop("comments")
        data_dict['pathways'] = [pathway.to_json_safe_dict() for pathway in self.pathways]
        data_dict['go_terms'] = {category: [term.to_json_safe_dict() for term in terms] for category, terms in self.go_terms.items()}
        data_dict['schema_version'] = CACHE_SCHEMA_VERSION
        return data_dict


class GenbankGeneFileComment(object):
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
            print('pathways')
            self._parse_pathways()

        if "GeneOntology" in self.headings:
            print('ontologies')
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