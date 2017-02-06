from bottle import route, run, template, static_file, view,\
                   post, request

def chargeandprop(aa_seq):
    """ Calculates protein net charge and charged AA proportion
    """
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

@route('/')
def index():
    return static_file('protchargeformbottle.html', root='views/')

@route('/css/<filename>')
def css_static(filename):
    return static_file(filename, root='css/')

@post('/protcharge')
@view('result')
def protcharge():
    seq = request.forms.get('aaseq', 'QWERTYYTREWQRTYEYTRQWE')
    prop = request.forms.get('prop','n')
    title = request.forms.get('title', 'No title')
    charge, propvalue = chargeandprop(seq)
    return {'seq': seq, 'prop': prop, 'title': title,
            'charge': round(charge, 3), 'propvalue': propvalue}

run(host='localhost', port=8000)
