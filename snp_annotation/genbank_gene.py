import re
import json

from bs4 import BeautifulSoup

HEADINGS = [
    "Pathways",
    "GeneOntology",
]

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
        self.biosytems = None
        self.go_terms = None

        if "xml" in opts:
            self._parse_xml(data)
        elif "json" in opts:
            self._from_json(data)

    def _parse_xml(self, data):
        self.parser = BeautifulSoup(data)
        self.gene_symbol = to_text_strip(self.parser.find("gene-ref_locus"))
        self.gene_ref_tag = to_text_strip(self.parser.find("gene-ref_locus-tag"))
        self.comments = [GenbankGeneFileComment(x.parent,self, xml=True) for x in self.parser.find_all('gene-commentary_heading') 
            if to_text_strip(x) in HEADINGS]


    def _from_json(self, json_dict):
        pass

    def _to_json_safe_dict(self):
        data_dict = self.__dict__
        data_dict.pop("parser")

def to_text(tag):
    return tag.get_text()

def to_text_strip(tag):
    return tag.get_text().strip()

def to_int(tag):
    return int(tag.get_text().strip())

def to_attr_value(tag):
    return tag.attrs['value']

class GenbankGeneFileComment(object):
    """docstring for GenbankGeneFileComment"""
    def __init__(self, data, owner, **opts):
        self.owner = owner
        self.opts = opts
        self.verbose = opts.get('verbose', False)

        self.substructures = []

        self.go_terms = {}

        if "xml" in opts:
            self._parse_xml(data)
        elif "json" in opts:
            self._from_json(data)

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
        subcomponents = [x for x in self.parser.find_all("gene-commentary")]
        subcomponents_text = [to_text_strip(x.find("gene-commentary_text")) for x in subcomponents]
        subcomponents_url  = [to_text_strip(x.find("other-source_url")) for x in subcomponents]
        print(subcomponents_text)
        print(subcomponents_url)
        self.contents.extend(subcomponents_text)
        self.urls.extend(subcomponents_url)


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
        
        self.go_terms.update(go_terms)




    def __repr__(self):
        rep = "Comment(Headings: %(headings)s Content: %(contents)s)"
        return rep


class Pathway(object):
    def __init__(self, name, url):
        self.name = name
        self.url = url

class GOTerm(object):
    def __init__(self, db, id, term):
        self.db = db
        self.id = id
        self.term = term

    def __repr__(self):
        rep = "GOTerm(%(db)s|%(id)s|%(term)s)" % self.__dict__
        return rep