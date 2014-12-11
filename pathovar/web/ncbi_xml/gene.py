import re
import json
from collections import defaultdict
from copy import copy

from pathovar.web.ncbi_xml import ET, to_text, to_text_strip, to_int, to_attr_value



HEADINGS = [
    "Pathways",
    "GeneOntology",
    "Interactions",
]

class GenBankGeneFile(object):
    CACHE_SCHEMA_VERSION = '0.3.8'
    def __init__(self, data, **opts):
        self.parser = None
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.gid = None
        self.gene_db_id = None
        self.gene_ref_tag = None
        self.gene_prot_name = None
        self.gene_symbol = None
        self.comments = None
        self.pathways = []
        self.go_terms = defaultdict(list)

        if "xml" in opts:
            self._parse_xml(data)
        elif "json" in opts:
            self._from_json(data)

    def _parse_xml(self, data):
        self.parser = ET.fromstring(data)
        self.gene_db_id = self.parser.find(".//Gene-track_geneid")
        if self.gene_db_id is not None:
            self.gene_db_id = self.gene_db_id.text
        self.gene_symbol = self.parser.find(".//Gene-ref_locus")
        if self.gene_symbol is not None:
            self.gene_symbol = to_text_strip(self.gene_symbol)
        self.gene_ref_tag = self.parser.find(".//Gene-ref_locus-tag")
        if self.gene_ref_tag is not None:
            self.gene_ref_tag = to_text_strip(self.gene_ref_tag)
        self.gene_prot_name = self.parser.find(".//Prot-ref_name_E")
        if self.gene_prot_name is not None:
            self.gene_prot_name = to_text_strip(self.gene_prot_name)
        self.comments = [GenBankGeneFileComment(x,self, xml=True) for x in self.parser.findall('.//Gene-commentary_heading/..') 
            if to_text_strip(x) in HEADINGS]
        for comment in self.comments:
            self.pathways.extend(comment.pathways)
            for key, value in comment.go_terms.items():
                self.go_terms[key].extend(value)

    def _from_json(self, json_dict):
        self.gene_db_id = json_dict['gene_db_id']
        self.gene_symbol = json_dict["gene_symbol"]
        self.gene_ref_tag =  json_dict["gene_ref_tag"]
        self.gene_prot_name = json_dict["gene_prot_name"]
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
        self.headings = [to_text_strip(h) for h in self.parser.findall(".//Gene-commentary_heading")]

        if "Pathways" in self.headings:
            self._parse_pathways()

        if "GeneOntology" in self.headings:
            self._parse_gene_ontologies()

    def _parse_pathways(self):
        pathways = []
        subcomponents = [x for x in self.parser.findall(".//Gene-commentary")]
        for subcomponent in subcomponents:
            text = to_text_strip(subcomponent.find(".//Gene-commentary_text"))
            db = to_text_strip(subcomponent.find(".//Dbtag_db"))
            id = to_text_strip(subcomponent.find(".//Object-id_str"))
            url = to_text_strip(subcomponent.find(".//Other-source_url"))
            pathway = Pathway(text, db, id, url)
            pathways.append(pathway)
        self.pathways = pathways

    def _parse_gene_ontologies(self):
        subcomponents = [x for x in self.parser.findall(".//Gene-commentary")]
        subcomponents_categories = [x.find(".//Gene-commentary_label") for x in subcomponents if x.find(".//Gene-commentary_label")]
        subcomponents_categories = {to_text_strip(x): x.parent for x in subcomponents_categories}

        go_terms = dict()
        for category_name in subcomponents_categories:
            category = subcomponents_categories[category_name]
            entries = category.findall(".//Gene-commentary")
            entry_terms = []
            for ent in entries:
                db = ent.find(".//Dbtag_db").text
                id = ent.find(".//Object-id_id").text
                term = ent.find(".//Other-source_anchor").text
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

class GenBankGeneToBioSystem(object):
    CACHE_SCHEMA_VERSION = "0.3"

    def __init__(self, data, **opts):
        self.parser = None
        self.opts = opts
        self.gene_id = None
        self.biosystem_ids = None

        if "xml" in self.opts:
            self._parse_xml(data)
        elif "json" in self.opts:
            self._from_json(data)

    def _parse_xml(self, data):
        self.parser = ET.fromstring(data)
        self.gene_id = self.opts['gene_id']
        self.biosystem_ids = list(set([link.text for link in self.parser.findall(".//Id") if link.text != self.gene_id]))


    def _from_json(self, json_dict):
        self.biosystem_ids = json_dict['biosystem_ids']
        self.gene_id = json_dict['gene_id']

    def to_json_safe_dict(self):
        data_dict = self.__dict__
        data_dict.pop('parser', None)
        data_dict['schema_version'] = GenBankGeneToBioSystem.CACHE_SCHEMA_VERSION
        return data_dict


class GenBankBioSystemFile(object):
    CACHE_SCHEMA_VERSION = "0.2b"
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
        self.parser = ET.fromstring(data)
        self.system_names = [to_text_strip(name) for name in self.parser.findall(".//System_names_E")]
        self.system_id = to_text_strip(self.parser.find(".//Sys-id_bsid"))
        self.system_description = to_text_strip(self.parser.find(".//System_description"))
        self.external_url = to_text_strip(self.parser.find(".//System_recordurl"))
        self.external_accession = to_text_strip(self.parser.find(".//System_externalaccn"))
        self.system_categories = [to_text_strip(cat) for cat  in self.parser.findall(".//System_category_E")]

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

if __name__ == '__main__':
    import sys
    file_name = sys.argv[1]
    data = ''.join(open(file_name).readlines())
    ggf = GenBankGeneFile(data, xml = True, verbose = True)
    import IPython
    print("GeneFile stored in local variable `ggf`")
    IPython.embed()
