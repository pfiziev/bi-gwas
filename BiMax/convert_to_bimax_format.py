import json
import os
import sys

__author__ = 'pf'


if __name__ == '__main__':

    input_fname = sys.argv[1]
    print input_fname

    data = json.load(open(input_fname))
    outf = open(os.path.split(input_fname)[1]+'.bimax', 'w')

    outf.write('%d %d %d %d\n' % (len(data['cases']), len(data['cases'][0]), 50, 50))

    outf.write('\n'.join(' '.join(str(0 if v == 0 else 1) for v in row) for row in data['cases']))

    outf.close()