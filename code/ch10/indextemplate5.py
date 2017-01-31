from bottle import route, run, template

@route('/greets/<name>')
def shows_greeting(name):
    return template('new', **{'name':name})

run(host='localhost', port=8000)
