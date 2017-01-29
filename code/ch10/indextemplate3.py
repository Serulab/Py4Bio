from bottle import route, run, template

@route('/greets/<username>')
def shows_greeting(username):
    if username[0].isalpha():
        msg = 'Hello {0}!'.format(username)
    else:
        msg = "Your username must can't start with a number"
    return template('main_template3', **{'msg':msg})

run(host='localhost', port=8000)
