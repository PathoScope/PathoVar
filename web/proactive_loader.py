import multiprocessing
import time

from pathovar.web import entrez_eutils
from pathovar.snp_annotation.locate_variant import GenBankFeatureFile

# Steals the functionality of EntrezAnnotationMapper without requiring a VCF to parse
# It is not safe to call any functions other than self.find_annotations_by_gene_id 
class AnnotationFetcher(locate_variant.EntrezAnnotationMapper):
    def __init__(self):
        self.entrez_handle = entrez_eutils.EntrezEUtilsDriver(verbose=False)
        self.annotation_cache = dict()
        self.verbose = False

def fetch_annotation_file(gi, pipe_to_manager):
    anno = AnnotationFetcher()
    data = anno.find_annotations_by_gene_id(gi, {"verbose":True})
    pipe_to_manager.send(data.to_json_safe_dict())


# Manage requests to EntrezEUtils to make sure that no more than 3 requests are made per second
# so that the IP address is not black-listed.
class EntrezRequestManager(object):
    READY = 0
    WAITING = 1
    
    def __init__(self, gis_to_request, max_processes = 5):
        self.gis_to_request = gis_to_request
        self.processes = []
        self.pipes = []
        self.results = []
        self.max_processes = max_processes
        self.last_request = time.time()
        self.state = READY

    def start_process(self, gi):
        if time.time() - self.last_request > 100: # Need to convert to seconds here. Look up the unit
            time.sleep(3)
        if len(self.processes) > self.max_processes:
            self.manage_processes()
            
        parent_conn, child_conn = multiprocessing.Pipe()
        proc = multiprocessing.Process(target=fetch_annotation_file, args=[gi,child_conn])
        proc.start()
        self.processes.append(proc)
        self.pipes.append(parent_conn)

    def manage_processes(self):
        pass

    def wait_for_process(self):
        for i,proc in enumerate(self.processes):
            if not proc.is_alive():
                self.results.append(self.pipes[i].recv())
                self.processes.pop(i)
                self.pipes.pop(i)


