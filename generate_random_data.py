__author__ = 'pf'
from pprint import pformat
import random

__author__ = 'pf'

SNPs = 20
CASES = 20
CONTROLS = 20
BICLUSTERS = 2
BI_MAX_SNPs = 5
BI_MAX_INDs = 5

def matrix(n, m):
    return [[0 for j in xrange(m)] for i in xrange(n)]

if __name__ == '__main__':
    case_m = matrix(SNPs, CASES)
    cont_m = matrix(SNPs, CONTROLS)
    for bc in xrange(BICLUSTERS):
        bc_i = random.randint(BI_MAX_INDs - 2, BI_MAX_INDs)
        bc_s = random.randint(BI_MAX_SNPs - 2, BI_MAX_SNPs)
        case_inds = random.sample(xrange(CASES), bc_i)
        case_snps = random.sample(xrange(SNPs), bc_s)
        for i in case_snps:
            for j in case_inds:
                case_m[i][j] = 1
    print pformat(case_m)
#       control_inds = random.sample(xrange(CONTROLS))

