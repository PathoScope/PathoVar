# Standard Library Modules
import json
import os
import re

TEMPLATE_PATH = 'template.html'

def embed_data(template_file, data_file, output_file = None):
    data = ''.join(open(data_file).readlines())
    template = open(template_file).readlines()
    if output_file is None:
        name, ext = os.path.splitext(template_file)
        name += '.emb' + ext
        output_file = name
    output_fh = open(output_file,'w')
    for line in template:
        if re.search(r"script\s*id='data-hook'", line):
            output_fh.write("<script async=true id='data-hook'>Organisms=" + data + "</script>")
        else:
            output_fh.write(line)
    output_fh.close()
    return output_file

def build_report(data_file, **opts):
    output_file = os.path.dirname(data_file) + os.sep + \
        os.path.splitext(os.path.basename(data_file))[0] + '_report.html'
    return embed_data(TEMPLATE_PATH, data_file, output_file)