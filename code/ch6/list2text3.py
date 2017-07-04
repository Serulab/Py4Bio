def commandline(name, **parameters):
    line = ''
    for item in parameters:
        line += ' -{0} {1}'.format(item, parameters[item])
    return name + line
