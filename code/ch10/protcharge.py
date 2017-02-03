#!/usr/bin/env python
import cgi, cgitb

def chargeandprop(aa_seq):
    protseq = aa_seq.upper()
    charge = -0.002
    cp = 0
    aa_charge = {'C':-.045,'D':-.999,'E':-.998,'H':.091,
                 'K':1,'R':1,'Y':-.001}
    for aa in protseq:
        charge += aa_charge.get(aa, 0)
        if aa in aa_charge:
            cp += 1
    prop = float(cp)/len(aa_seq)*100
    return (charge, prop)

cgitb.enable()
print('Content-Type: text/html\n')
form = cgi.FieldStorage()
uname = form.getvalue('username','NN')
seq = form.getvalue('seq', 'QWERTYYTREWQRTYEYTRQWE')
prop = form.getvalue('prop', 'n')
jobtitle = form.getvalue('title','No title')
charge, propvalue = chargeandprop(seq)
print('<html><body>Job title:{0}<br/>'.format(jobtitle))
print('Your sequence is:<br/>{0}<br/>'.format(seq))
print('Net charge: {0}<br/>'.format(charge))
if prop == 'y':
 print('Proportion of charged AA: {0:.2f}<br/>'
       .format(propvalue))
print('</body></html>')
