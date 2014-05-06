import re
import json
from time import sleep, time
from copy import copy
from collections import OrderedDict, defaultdict

from pathovar.utils.fasta_utils import SequenceRecord
from pathovar.web import entrez_eutils
from pathovar.web.ncbi_xml import ET, to_text, to_text_strip, to_int, to_attr_value, try_tag, take_def

## GenBankFeatureFile
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the XML definition of a GenBank flat file.
class GenBankFeatureFile(object):
    CACHE_SCHEMA_VERSION = '0.3.10'
    def __init__(self, data, **opts):
        self.opts = opts
        self.verbose = opts.get('verbose', False)
        self.mol_type = opts.get('mol_type', 'nucl')
        timer = time()
        self.gid = None
        self.accession = None
        self.parser = None
        self.chromosome = None
        self.org_name = None
        self.genetic_code = None
        self.entries = {}
        self.genome_entry = None
        self.features = []
        self.sorted_genes = []
        self.last_entry_ind = 0
        self.db_tag_data = dict()
        if 'xml' in opts:
            self._parse_xml(data)
        elif 'json' in opts:
            self._from_json(data)

    def _parse_xml(self, xml):
        timer = time()
        self.parser = ET.fromstring(xml)
        if self.verbose: print("XML Digested (%s sec)" % str(time() - timer))
        if self.verbose: print("Searching for Chromosome")
        self.gid = self.parser.find(".//Seq-id_gi").text
        if self.mol_type == 'nucl':
            self.org_name = self.parser.find(".//Org-ref_taxname").text
            org_mod_name = self.parser.find(".//OrgName_mod")
            if org_mod_name is not None:
                org_mod_subtype = org_mod_name.find(".//OrgMod_subtype")
                org_mod_subtype = to_attr_value(org_mod_subtype, 'value') if org_mod_subtype is not None else ""

                org_mod_subname = org_mod_name.find(".//OrgMod_subname")
                org_mod_subname = org_mod_subname.text if org_mod_subname is not None else ""
                self.org_name += " " + org_mod_subtype + " " + org_mod_subname
            
            subsource = self.parser.find(".//BioSource_subtype")
            if subsource is not None:
                subtype = subsource.findall(".//SubSource_subtype")
                subtype_name = subsource.findall(".//SubSource_name")
                #subtype = to_attr_value(subtype, 'value') if subtype is not None else ''
                #subtype_name = subtype_name.text if subtype_name is not None else ''
                for i in range(len(subtype_name)):
                    #subtype_i = subtype[i]
                    subtype_name_i = subtype_name[i]
                    #subtype_i = to_attr_value(subtype_i, 'value') if subtype_i is not None else ''
                    subtype_name_i = subtype_name_i.text if subtype_name_i is not None else ''
                    self.org_name += ' ' + subtype_name_i
            db_tags = self.parser.find(".//Org-ref_db")
            if db_tags is not None:
                db_tags = db_tags.findall(".//Dbtag")
                self.db_tag_data = dict()
                for tag in db_tags:
                    db_name = tag.find(".//Dbtag_db")
                    db_name = db_name.text if db_name is not None else ""
                    tag_data = to_text_strip(tag.find(".//Dbtag_tag"))
                    self.db_tag_data[db_name] = tag_data


            self.genetic_code = int(self.parser.find(".//OrgName_gcode").text)
            self.accession = self.parser.find(".//Textseq-id_accession").text + '.' + self.parser.find(".//Textseq-id_version").text
        #if self.mol_type == 'nucl':
            self.chromosome = self.parser.findall('.//IUPACna')
            if self.chromosome:
                self.chromosome = ''.join(map(to_text_strip, self.chromosome))
                if self.verbose: print("Found, %d bp" % len(self.chromosome))
            chromosome_len_parity = sum(map(to_int, self.parser.findall('.//Seq-literal_length')))

            if len(self.chromosome) != chromosome_len_parity and chromosome_len_parity != 0:
                if self.verbose: print("Chromosome Length Parity Error (%d != %d). Sequence Data Missing. Attempting to fix" % (len(self.chromosome), chromosome_len_parity))
                eutils_handle = entrez_eutils.EntrezEUtilsDriver(**self.opts)
                data = eutils_handle.find_nucleotides_by_gene_id(self.gid)
                fasta_seq = (data.split('\n'))
                replace_chromosome = ''.join(fasta_seq[1:])
                if len(replace_chromosome) == chromosome_len_parity:
                    if self.verbose: print("Chromosome Length Parity Fixed")
                    self.chromosome = replace_chromosome
                else:
                    raise GenBankFileChromosomeLengthMismatch("Could not resolve Chromosome Length Parity Error")

        if self.verbose: print("Gathering Entries and Features")
        self.entries = {ent.gid : ent for ent in map(lambda x: GenBankSeqEntry(x, self, xml = True), self.parser.findall(".//Seq-entry")) }
        self.genome_entry = self.entries[self.gid]
        # Create Features for each Seq-feat tag, both in the Genome level sequence and the individuals
        self.features = map(lambda x: GenBankFeature(x, self), self.parser.findall(".//Seq-feat"))

        # Features are subsets of Entries. It would be a good idea to compress them to a single entity
        # later. Entries capture finer resolution details about a particular gene
        for feature in self.features:
             if feature.gid in self.entries:
                 entry = self.entries[feature.gid]
                 feature.title = entry.title
                 entry.strand = feature.strand
        
        genome_level_features = {}
        sequence_level_features = []
        
        print("Mapping Features Over Genomic Position")
        for feat in self.features:
            if "complete genome" in feat.title:
                if feat.start not in genome_level_features:
                    genome_level_features[feat.start] = feat
                else: 
                    genome_level_features[feat.start].merge_features(feat)
            else:
                sequence_level_features.append(feat)
        
        # Each feature occurs multiple times in the file, redundant with its multiple regions. The complete
        # genomic span is the largest span. This works for single-span entities. 
        
        if self.verbose: print("Computing Genomic Coordinates")
        feature_dict = {}
        for feat in sequence_level_features:
            # Discard feature if its gid matches the genome gid
            if feat.gid == self.gid: continue
            if feat.starts:
                if feat.gid not in feature_dict:
                    feature_dict[feat.gid] = feat
                else:
                    if (feat.end - feat.start) > (feature_dict[feat.gid].end - feature_dict[feat.gid].start):
                        feature_dict[feat.gid] = feat
                    if self.mol_type == 'nucl' and feat.end > len(self.chromosome):
                        raise GenBankFileChromosomeLengthMismatch("%r exceeds chromosome size" % feat)
        
        sequence_level_features = feature_dict.values()
        by_start = {f.start: f for f in sequence_level_features}
        
        final_features = {}
        # Reduce all features to a non-redundant set based on
        # start position and identity. 
        for pos in genome_level_features:
            try:
                final_features[by_start[pos].gid] = genome_level_features[pos]
                final_features[by_start[pos].gid].gid = by_start[pos].gid
                genome_level_features[pos].gid = by_start[pos].gid
            except KeyError, e:
                # Using a compound of title + locus tag in place of absent gid
                genome_level_features[pos].title = genome_level_features[pos].gene_ref_tag
                final_features[genome_level_features[pos].title] = genome_level_features[pos]
                final_features[genome_level_features[pos].title].gid = genome_level_features[pos].title


        merged_features = final_features.values()
        # If there were no genome level features labeled, then nothing is kept.
        # Instead, keep all of the merged 'sequence level features'
        if len(merged_features) == 0:
            merged_features = sequence_level_features
        self.features = merged_features

        # Using the genomic position data just computed, update coordinate information for the 
        # related entries
        for feat in self.features:
            try:
                entry = self.entries[feat.gid]
                entry.update_genome_position(feat)
            except:
                if feat.gid is not None or feat.gene_ref_tag is not None:
                    entry = feat.upgrade_to_entry()
                    self.entries[entry.gid] = entry
        
        # Remove whole-reference entry
        self.entries.pop(self.gid)
        self.sorted_genes = sorted([gene for gid, gene in self.entries.items()], key=lambda x: x.start)

    def _from_json(self, json_dict):
        timer = time()
        self.gid = json_dict['gid']
        self.accession = json_dict['accession']
        self.org_name = json_dict['name']
        self.chromosome = json_dict['chromosome']
        self.genetic_code = json_dict.get('genetic_code', None)
        self.entries = {k:GenBankSeqEntry(v, self, **self.opts) for k,v in json_dict['entries'].items()}
        self.sorted_genes = json_dict["sorted_genes"]
        self.sorted_genes = [self.entries[gid] for gid in self.sorted_genes]
        self.db_tag_data = json_dict['db_tag_data']
        if self.verbose: print("Loading from JSON Complete (%rs)" % (time() - timer))


    def to_json_safe_dict(self):
        data_dict = {}
        data_dict['gid'] = self.gid
        data_dict['accession'] = self.accession
        data_dict['chromosome'] = self.chromosome
        data_dict['schema_version'] = GenBankFeatureFile.CACHE_SCHEMA_VERSION
        data_dict["name"] = self.org_name
        data_dict['genetic_code'] = self.genetic_code
        if self.entries is None: 
            print("Entries is None")
        data_dict['entries'] = {
            k: v.to_json_safe_dict() for k,v in 
            self.entries.items() if v.title or "complete genome" not in v.title
            }
        bang = False
        for x in self.sorted_genes:
            if x.gid is None:
                print(x.__dict__)
                bang = True
        if bang: exit(1)
        data_dict['sorted_genes'] = [x.gid for x in self.sorted_genes]
        data_dict['db_tag_data'] = self.db_tag_data
        return(data_dict)

    ##
    # 
    def locate_snp_site(self, snp, reset=True, forward_check = False):
        results = []
        if reset:
            self.last_entry_ind = 0
        for ind, entry in enumerate(self.sorted_genes[self.last_entry_ind:]):
            #if self.verbose: print(entry.start, entry.end)
            #if self.verbose: print("check %d >= %d and %d <= %d" % (snp_loc, entry.start, snp_loc, entry.end))
            if (snp.start >= entry.start and snp.start <= entry.end) or (snp.end >= entry.end and snp.start <= entry.end) \
            or (snp.start <= entry.start and snp.end >= entry.start):
                #if self.verbose: print("SNP Location %r mapped within %r" % ([snp.start, snp.end], entry))
                self.last_entry_ind += ind
                if(forward_check):
                    results.append(entry)
                    continue
                return entry
            if snp.start < entry.start and snp.end < entry.start:
                # return intergenic relative to this entry and the previous one
                upstream_id = None
                downstream_id = entry.gid
                if(self.last_entry_ind + ind != 0):
                    upstream_id = self.sorted_genes[ind - 1].gid
                entry = IntergenicEntry(upstream_id, downstream_id)
                if(forward_check):
                    results.append(entry)
                    continue
                return entry
            if snp.start < entry.start and snp.end > entry.end and forward_check:
                print('big step')
                results.append(entry)
        # Forward-Check mode expects a list result
        if forward_check:
            if len(results) > 0:
                return results
            else:
                return [IntergenicEntry(self.sorted_genes[-1].gid)]
        return IntergenicEntry(self.sorted_genes[-1].gid)

    def __repr__(self):
        rep = "GenBankFile(" + self.mol_type + "|" + ', '.join(map(repr, self.entries)) + ")"
        return rep

## GenBankFeature
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the substructure of a GenBank flat file related to a single sequence feature
class GenBankFeature(object):
    def __init__(self, parser, owner):
        self.parser = parser
        self.owner = owner
        self.gid = to_text_strip(parser.find('.//Seq-id_gi'))
        self.components = None
        # Some features, especially in complex organisms, will have multi-part features. 
        # Capture all of that data
        self.starts = map(to_int, self.parser.findall('.//Seq-interval_from'))
        self.ends = map(to_int, self.parser.findall('.//Seq-interval_to'))

        self.title = None

        # and set the extrema to the global start and stop
        self.start = None
        self.end   = None

        if(self.starts):
            self.start = min(self.starts)
            self.end = max(self.ends)
        else:
            pass

        self.gene_symbol = self.parser.find(".//Gene-ref_locus")
        if self.gene_symbol is not None:
            self.gene_symbol = to_text_strip(self.gene_symbol)
        self.gene_ref_tag = self.parser.find(".//Gene-ref_locus-tag")
        if self.gene_ref_tag is not None:
            self.gene_ref_tag = to_text_strip(self.gene_ref_tag)
        self.gene_prot_name = self.parser.find(".//Prot-ref_name_E")

        self.strand = self.parser.find('.//Na-strand')
        if self.strand is not None:
            self.strand = self.strand.get('value')
            assert self.strand != ''
        else:
            self.strand = "?"

        self.is_rna = self.parser.find(".//RNA-ref")
        self.rna_type = None
        self.rna_desc = None
        if self.is_rna is not None:
            self.rna_type = to_attr_value(self.is_rna.find(".//RNA-ref_type"), 'value')
            self.rna_desc = try_tag(self.parser.find(".//RNA-ref_ext_name"), to_text)
            if self.rna_desc is None:
                self.rna_desc = try_tag(self.parser.find(".//Trna-ext_aa_ncbieaa"), to_text)
            if self.rna_desc is None:
                self.rna_desc = "miscRNA"

        self.is_rna = True if self.is_rna is not None else False

        #self.__dict__.pop('parser')
    
    def merge_features(self, other):
        merge_dict = dict()
        hold = self.components
        for key in self.__dict__:
           merge_dict[key] = take_def(key, self.__dict__, other.__dict__)
        self.__dict__ = merge_dict
        self.components = (hold, self, other)
        


    def upgrade_to_entry(self):
        upgrade_dict = dict()
        upgrade_dict["start"] = self.start
        upgrade_dict["end"] = self.end
        upgrade_dict["starts"] = self.starts
        upgrade_dict["ends"] = self.ends
        upgrade_dict["strand"] = self.strand
        upgrade_dict["gene_symbol"] = self.gene_symbol
        upgrade_dict["gene_ref_tag"] = self.gene_ref_tag
        upgrade_dict["gene_prot_name"] = self.gene_prot_name
        upgrade_dict["gid"] = self.gid
        upgrade_dict["accession"] = self.gid
        upgrade_dict["title"] = self.rna_desc if self.is_rna else ""
        upgrade_dict["annotations"] = {}
        upgrade_dict["comments"] = []
        upgrade_dict["amino_acid_sequence"] = None
        upgrade_dict["nucleotide_sequence"] = self.owner.chromosome[self.start:self.end] if (self.owner.mol_type == 'nucl') else None
        upgrade_dict["is_rna"] = self.is_rna
        return GenBankSeqEntry(upgrade_dict, self.owner, json=True)


    def __repr__(self):
        rep = "GenBankFeature(gi|%(gid)s, %(start)r, %(end)r, %(title)s)" % self.__dict__
        return rep

## GenBankSeqEntry
# XML structure parser and annotation extraction object. Uses BeautifulSoup to parse
# the substructure of a GenBank flat file related to a single sequence entry
class GenBankSeqEntry(object):
    def __init__(self, data, owner, **opts):
        self.parser = None
        self.owner = owner

        self.gid = None
        self.accession = None
        
        self.gene_ref_tag = None
        self.gene_symbol = None
        self.protein_name = None

        # The coordinates in the seq-entry itself are relative to their own start and stop points. 
        # The genomic coordinates must be inferred from the genomic seq-feat tags
        self.starts = None
        self.ends   = None
        self.start  = None
        self.end    = None
        
        self.strand = None
        
        self.is_partial = None
        self.is_pseudo = None
        self.is_rna = False
        self.is_intergenic = False

        self.amino_acid_sequence = None
        # Prune later once starts and ends are set
        self.nucleotide_sequence = None 

        self.title = None
        self.annotations = {}
        self.comments = []
        
        self._id = None
        self._desc = None
        self._inst = None
        self._annot = None

        if 'xml' in opts:
            self._parse_xml(data)
        elif 'json' in opts:
            self._from_json(data)

    def _parse_xml(self, parser):
        self.parser = parser

        self.gid = to_text_strip(parser.find('.//Seq-id_gi'))

        self.accession = to_text_strip(parser.find('.//Textseq-id_accession'))
        
        self.strand = self.parser.find('.//Na-strand')
        if self.strand is not None:
            self.strand = self.strand.get('value')
            assert self.strand != ''
        else:
            self.strand = "?"

        self.is_partial = self.parser.find(".//Seq-feat_partial")
        if self.is_partial is not None: 
            self.is_partial = self.is_partial.get('value')
        else:
            self.is_partial = False

        self.is_pseudo = self.parser.find(".//Seq-feat_pseudo")
        if self.is_pseudo is not None: 
            self.is_pseudo = self.is_pseudo.get('value')
        else:
            self.is_pseudo = False

        self.amino_acid_sequence = map(to_text, parser.findall('.//IUPACaa'))
        # Prune later once starts and ends are set
        self.nucleotide_sequence = self.owner.chromosome

        self.title = to_text_strip(parser.find('.//Seqdesc_title'))
        #self.annotations = {ann.name: ann for ann in map(lambda d: GenBankAnnotation(d, xml=True, verbose=self.owner.verbose), parser.findall("seq-annot"))}
        self.annotations = {ann.name: ann for ann in map(lambda d: GenBankAnnotation(d, xml=True, verbose=self.owner.verbose), parser.findall(".//Seq-feat"))}

        self.comments = map(to_text_strip, parser.findall('.//Seq-feat_comment'))


    def _from_json(self, json_dict):
        try:
            self.gid = json_dict['gid']
            self.accession = json_dict['accession']
            
            self.strand = json_dict['strand']
            self.amino_acid_sequence = json_dict['amino_acid_sequence']
            self.nucleotide_sequence = json_dict['nucleotide_sequence']
            
            self.title = json_dict['title']
            self.comments = json_dict['comments']
            
            self.start =  json_dict["start"]
            self.end = json_dict["end"]
            self.starts = json_dict["starts"]
            self.ends =  json_dict["ends"]
            
            self.is_partial = json_dict.get('is_partial', False)
            self.is_pseudo = json_dict.get('is_pseudo', False)

            self.gene_symbol = json_dict["gene_symbol"]
            self.gene_ref_tag = json_dict["gene_ref_tag"]
            self.gene_prot_name = json_dict.get("gene_prot_name", None)

            # Only ever set through an upgrade from Feature, which is passed by dictionary
            self.is_rna = json_dict["is_rna"]

            self.annotations = {k:GenBankAnnotation(v, json=True) for k,v in json_dict['annotations'].items()}
        except KeyError, e:
            print("Error finding key %s for %s")

    def update_genome_position(self, feature):
        self.starts = feature.starts
        self.ends   = feature.ends  
        self.start  = feature.start 
        self.end    = feature.end
        self.gene_symbol = feature.gene_symbol
        self.gene_ref_tag = feature.gene_ref_tag
        self.gene_prot_name = feature.gene_prot_name

        if self.owner.mol_type == 'nucl':
            self.nucleotide_sequence = self.owner.chromosome[self.start:self.end]

    def __repr__(self):
        rep = "GenBankSeqEntry(gi=%(gid)s, Title=%(title)s, ACC=%(accession)s)" 
        return rep % self.__dict__

    def to_info_field(self):
        rep = "(gi:%(gid)s|title:%(title)s|acc:%(accession)s|strand:%(strand)s" % self.__dict__
        annos = map(lambda x: x.to_info_field(), self.annotations.values())
        rep_annos = '||Annotations|' + '|'.join(annos)
        rep_comments = '||Comments|' + '|'.join(self.comments)
        if re.search(r'resistance', rep + rep_annos + rep_comments):
            rep += '||RESISTANCE'
        rep += ')'
        rep = re.sub(r'\s', '_', rep)
        return rep

    def to_json_safe_dict(self):
        data_dict = copy(self.__dict__)
        # Remove fields that don't serialize.
        unsafe_fields = ["parser", "_id","_desc","_inst","_annot","owner"]
        for f in unsafe_fields:
            try:
                data_dict.pop(f)
            except KeyError:
                pass
            # Translate nested objects to dictionaries
        data_dict["annotations"] = {
                k:v.__dict__ for k,v in 
                data_dict["annotations"].items()
            }


        return data_dict

class IntergenicEntry(object):
    def __init__(self, upstream = None, downstream = None):
        self.upstream_id = upstream
        self.downstream_id = downstream

    @property
    def gid(self):
        return [self.upstream_id, self.downstream_id]

    def is_intergenic(self):
        return True

    def to_info_field(self):
        return "(Intergenic:%s~%s)" % (self.upstream_id, self.downstream_id)

class AnnotationExtension(OrderedDict):
    def __init__(self,*args, **kwargs):
        OrderedDict.__init__(self, *args, **kwargs)

    def __eq__(self, other_dict):
        for k in self.keys():
            if self[k] != other_dict[k]: return False
        print(self)
        print(other_dict)
        raw_input("Next")
        return True

    def __cmp__(self, other_dict):
        for k in self.keys():
            if self[k] != other_dict[k]: return False
        print(self)
        print(other_dict)
        raw_input("Next")
        return True

    def __repr__(self):
        rep = ''
        ext_type = self.get("ext_type", 'anno_ext')
        rep += ext_type + '('
        vals = []
        for key in self:
            if key == "ext_type": continue
            vals.append(key + ': ' + str(self[key]))
        rep += ', '.join(vals) + ')'
        return rep

    def __str__(self):
        return repr(self)

class GenBankAnnotation(object):
    def __init__(self, data, **opts):
        self.id = None
        self.start = []
        self.end = []
        self.name = None
        self.comments = []
        self.regions = []
        if 'xml' in opts:
            self._parse_xml(data)
        elif 'json' in opts:
            self._from_json(data)

    def _parse_xml(self, data):
        self.id = map(to_text_strip, data.findall(".//Seq-feat_dbxref"))
        self.start = map(lambda x: int(x.text.replace(' ','')), data.findall(".//Seq-interval_from"))
        self.end = map(lambda x: int(x.text.replace(' ','')), data.findall(".//Seq-interval_to"))
        
        self.regions = []#list(map(lambda x: x.text.strip(), data.findall("seqfeatdata_region")))
        self.name = ', '.join(map(to_text_strip, data.findall(".//Gene-ref") + data.findall(".//Prot-ref_name_e")))
        self.comments = (map(to_text_strip, data.findall(".//Seq-feat_comment")))
        if self.name == '' and len(self.comments) > 0:
            self.name = self.comments[0]
        raw_extensions = data.findall(".//Seq-feat_ext")
        for raw_ext in raw_extensions:
            obj_type = to_text_strip(raw_ext.find(".//Object-id_str"))
            obj_data = map(to_text_strip, raw_ext.findall(".//User-object_data")[0].findall(".//User-field_data"))
            ext = AnnotationExtension()
            ext = self.process_extension(raw_ext)
            ext['ext_type'] = obj_type
            self.regions.append(ext)

    def _from_json(self, json_dict):
        self.id = json_dict['id']
        self.start = json_dict['start']
        self.end = json_dict['end']
        self.name = json_dict['name']
        self.regions = json_dict['regions']
        self.comments = json_dict['comments']

    def __repr__(self):
        rep = "GenBankAnnotation(%(id)s|Starts=%(start)s, Ends=%(end)s, name:%(name)s, regions:%(regions)s)" % self.__dict__
        return rep

    def process_extension(self, raw_ext):
        labels = map(to_text_strip, raw_ext.findall(".//User-field_label"))
        data = raw_ext.findall(".//User-field_data")
        values = []
        for datum in data:
            text = to_text_strip(datum)
            if len(text) == 0:
                text = str(datum[0].get('value'))
                if text == '[]':
                    text = ''
            values.append(text)
        ext = AnnotationExtension()
        for i in xrange(len(labels)):
            ext[labels[i]] = values[i]

        return ext

    def add_regions(self, other_regions):
        for other_region in other_regions:
            if all([region != other_region for region in self.regions]):
                self.regions.append(other_region)

    def to_info_field(self):
        rep = '(starts:%(start)s|ends:%(end)s|' % self.__dict__
        for region in self.regions:
            buff = []
            region = copy(region)
            ext_type = region.pop('ext_type')
            pos_from = region.pop('from', None)
            pos_to   = region.pop('to', None)
            for key, value in region.items():
                buff.append(key + ":" + str(value))
            # prepend indices
            if pos_from:
                buff = ['from:' + str(pos_from).strip(), 'to:' + str(pos_to).strip()] + buff
            buff = '|'.join(buff)
            buff = ext_type + '(' + buff
            rep += buff
            rep += ')'
        rep += ')'
        return rep


class UnknownAnnotationException(Exception):
    pass
class GenBankFileChromosomeLengthMismatch(Exception):
    pass

if __name__ == '__main__':
    import sys, IPython
    file_name = sys.argv[1]
    data = ''.join(open(file_name).readlines())
    print("Parsing %s" % file_name)
    gbf = GenBankFeatureFile(data, mol_type = 'nucl' , xml = True, verbose = True)
    print("Data stored in local variable `gbf`")
    IPython.embed()



