import sys
try:
    x = open('random_filename')
except:
    a, b = sys.exc_info()[:2]
    print('Error name: {}'.format(a.__name__))
    print('Error code: {}'.format(b.args[0]))
    print('Error message: {}'.format(b.args[1]))
