from bottle import route, run

@route('/')
def index():
    """Display home page"""
    return '<h2>Hello World!</h2>'

run(host='localhost', port=8000)
