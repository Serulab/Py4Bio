#!/usr/bin/env python
import cgi
print("Content-Type: text/html\n")
form = cgi.FieldStorage()
name = form.getvalue("username","NN")[:10]
print("<html><head><title>A CGI script</title></head>")
print("<body><h2>Hello {0}</h2></body></html>".format(name))

