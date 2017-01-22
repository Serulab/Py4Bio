#!/usr/bin/env python

import os
import subprocess
from tempfile import mkstemp

from bottle import route, post, run, template, static_file, request
from bottle import debug

debug(True)

@route('/')
def index():
    return static_file('index.tpl', root='views/')

@post('/muscle_result')
def result():
    iterations = request.forms.get('iterat','1')
    output_type = request.forms.get('output','FASTA')
    order = request.forms.get('outorder','group')
    sequence = request.forms.get('seq','')
    if not sequence:
        # Since the textarea is empty, check the uploaded file
        sequence = request.files.get('upfile').file.read()
        #import pdb;pdb.set_trace()
    badrequest = ''
    # Verify that the user entered valid information.
    if iterations not in set(('1', '4', '8', '10', '12',
                              '14', '16')):
        badrequest = 'iterations'
    valid_output = set(('html','fasta','msf','clw','clwstrict'))
    if output_type not in valid_output:
        badrequest = 'output'
    if order not in set(('group', 'stable')):
        badrequest = 'outorder'
    result_output = ''
    # Make a random filename for user entered data
    fi_name = mkstemp('.txt','userdata_')[1]
    with open(fi_name,'wb') as fi_fh:
        fi_fh.write(sequence)
    fo_name = mkstemp('.txt','outfile_')[1]
    erfh = open('err.log','w')
    cmd = ['muscle3.8.31_i86linux64', '-in', fi_name,
           '-out', fo_name, '-quiet', '-maxiters',
           iterations, '-{}'.format(output_type),
           '-{}'.format(order)]
    p = subprocess.Popen(cmd, stderr=erfh)
    # Wait until finished
    p.communicate()
    erfh.close()
    with open(fo_name) as fout_fh:
        result_output = fout_fh.read()
    if output_type != 'html':
        result_output = '<pre>{}</pre>'.format(result_output)
    d = {'bad_option':badrequest, 'result_output':result_output}
    return template('result', d)

@route('/css/<filename>')
def css_static(filename):
    return static_file(filename, root='css/')

@route('/js/<filename>')
def js_static(filename):
    return static_file(filename, root='js/')

run(host='localhost', port=8080)
