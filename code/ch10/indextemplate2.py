from bottle import route, run, template

@route('/greets/<username>')
def shows_greeting(username):
    return template('main_template2', **{'name':username})

run(host='localhost', port=8000)
