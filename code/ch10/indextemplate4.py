from bottle import route, run, template

@route('/forloop')
def forloop():
    items = list(range(10))
    return template('index4', **{'items':items})

run(host='localhost', port=8000)
