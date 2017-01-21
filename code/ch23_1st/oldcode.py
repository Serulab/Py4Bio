#!/usr/bin/env python

import cgi
import cgitb
import subprocess
import sys
import os
from tempfile import mkstemp

# Uncomment the following line when debugging
#cgitb.enable()

def badrequest(bad):
    """ Display an error message """
    print("<h1>Bad Request</h1>\n")
    print("Use the options provided in the form: %s"%bad)
    print("</body></html>")
    # Get out of here:
    return sys.exit()

print("Content-Type: text/html\n")

form = cgi.FieldStorage()
iterat = form.getvalue("iterat","4")
output = form.getvalue("output","html")
outorder = form.getvalue("outorder","group")
# Get sequence data from text area
seqs = form.getvalue("seq")
if not seqs:
    # Since the textarea is empty, check the uploaded file
    seqs = form.getvalue("upfile")

# Verify that the user entered valid information.
if iterat not in set(('1','4','8','10','12','14','16')):
    badrequest(iterat)
valid_output = set(('html','fasta','msf','clw','clwstrict'))
if output not in valid_output:
    badrequest(output)
if outorder not in set(('group', 'stable')):
    badrequest(outorder)

print "<html><head><title>A CGI script</title></head><body>"

# Make a random filename for user entered data
fi_name=mkstemp('.txt','userdata_',"/var/www/muscleweb/")[1]
# Open this random filename
fi_fh = open(fi_name,'w')
# Write the user entered sequences into this temporary file
fi_fh.write(seqs)
fi_fh.close()

# Make a random filename for user entered data
fo_name=mkstemp('.txt','outfile_',"/var/www/muscleweb/")[1]
erfh = open('err.log','w')
cmd = ['./muscle', '-in', fi_name, '-out', fo_name,
       '-quiet', '-maxiters', iterat, '-%s'%output,
        '-%s'%outorder]
# Uncomment to check the generated command
#print ' '.join(cmd)
# Run the program with user provided parameters
p = subprocess.Popen(cmd, stderr=erfh, cwd='/var/www/muscleweb')
# Wait until finished
p.communicate()  # Same result as os.waitpid(p.pid,0)
erfh.close()
# Remove the input file since the it is not needed anymore.
os.remove(fi_name)
fout_fh = open(fo_name)
if output=='html':
    print(fout_fh.read())
else:
    print('<pre>%s</pre>'%fout_fh.read())
fout_fh.close()
# Remove the output file
os.remove(fo_name)
print("</body></html>")
This code is part of the book "Python for Bioinformatics", by Sebastian Bassi (sbassi@genesdigitales.com). Return to home page.
