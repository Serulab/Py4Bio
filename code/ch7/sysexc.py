import sys
try:
    0/0
except:
    a,b,c = sys.exc_info()
    print('Error name: {0}'.format(a.__name__))
    print('Message: {0}'.format(b))
    print('Error in line: {}'.format(c.tb_lineno))
