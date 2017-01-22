#!/usr/bin/env python

import os
import subprocess
from tempfile import mkstemp

from bottle import route, post, run, static_file, request, view
from bottle import debug

debug(True)

@route('/')
def index():
    return static_file('index.tpl', root='views/')

@post('/muscle_result')
@view('result')
def result():
    iterations = request.forms.get('iterat','1')
    output_type = request.forms.get('output','FASTA')
    order = request.forms.get('outorder','group')
    sequence = request.forms.get('seq','')
    if not sequence:
        # If the textarea is empty, check the uploaded file
        sequence = request.files.get('upfile').file.read()
    badreq = ''
    # Verify that the user entered valid information.
    try:
        int(iterations)
    except ValueError:
        badreq = 'iterations'
    valid_output = ('html', 'fasta', 'msf', 'clw', 'clwstrict')
    if output_type not in valid_output:
        badreq = 'output'
    if order not in ('group', 'stable'):
        badreq = 'outorder'
    result_out = ''
    # Make a random filename for user entered data
    fi_name = mkstemp('.txt','userdata_')[1]
    with open(fi_name,'wb') as fi_fh:
        fi_fh.write(sequence)
    fo_name = mkstemp('.txt','outfile_')[1]
    with open('muscle3_error.log','w') as erfh:
        cmd = ['muscle3.8.31_i86linux64', '-in', fi_name,
               '-out', fo_name, '-quiet', '-maxiters',
               iterations, '-{}'.format(output_type),
               '-{}'.format(order)]
        p = subprocess.Popen(cmd, stderr=erfh)
        p.communicate()
    # Remove the intput file
    os.remove(fi_name)
    with open(fo_name) as fout_fh:
        result_out = fout_fh.read()
    if output_type != 'html':
        result_out = '<pre>{}</pre>'.format(result_out)
    # Remove the output file
    os.remove(fo_name)
    return {'bad_option':badreq, 'result_output':result_out}

@route('/css/<filename>')
def css_static(filename):
    return static_file(filename, root='css/')

run(host='localhost', port=8080)
