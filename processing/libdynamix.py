import re

def read_param(paramname, filename):
    '''
    reads a parameter (float value) from a file in the format
    parametername=parametervalue
    '''
    # bad code! bad! should be able to import in other scripts
    f = open(filename)
    lines = f.readlines()
    for line in lines:
        if line[0:len(paramname)] == paramname:
            #print "found %s" % paramname
            #print line
            param = float(re.split('[= \t\n]',line.strip())[1])
            f.close()
            return param
    f.close()
    print "ERROR [read_param]: didn't find %s" % paramname
    return 0

def file_len(filename):
    with open(filename, 'r') as f:
        for i,l in enumerate(f):
            pass
    return i + 1
