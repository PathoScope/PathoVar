
# Standard Library Modules
import json
import os
import re

# Internal Modules
from pathovar.utils import vcf_utils

##
# Builds a all-in-one-file report for the data in a 
class StaticReport(object):
    TEMPLATE_PATH = "visualize/static/template.html"
    RESOURCE_PATH = "visualize/static/resources/"

    link_re = re.compile(r'''<link [^\n]* href=["']([^"']+)["']>''')
    scrip_re = re.compile(r'''<script (id=["'](P<id>[^"']+)["'])? .*>''')

    def __init__(self, vcf_path, **opts):
        self.vcf_path = vcf_path
        self.opts = opts
        self.verbose = opts.get("verbose", False)
        self.template_buffer = None

    def slurp_template(self):
        self.template_buffer = open(self.TEMPLATE_PATH).readlines()
        



