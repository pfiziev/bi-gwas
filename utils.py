import sys
import datetime

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)


global_stime = datetime.datetime.now()
def elapsed(msg = None):
    print "[%s]" % msg if msg is not None else "+", "Last:" , datetime.datetime.now() - elapsed.stime, '\tTotal:', datetime.datetime.now() - global_stime

    elapsed.stime = datetime.datetime.now()

elapsed.stime = datetime.datetime.now()

def matrix(n, m):
    return [[0 for j in xrange(m)] for i in xrange(n)]
