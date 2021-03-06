"""
This script is RETIRED!
"""


import sys
import datetime
import math

__author__ = 'pf'
from pprint import pformat
import random



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


def bit_encode(mat):
    res = [0]*len(mat)
    for i, row in enumerate(mat):
        mask = 1
        for cell in row:
            if cell:
                res[i] += mask
            mask <<= 1

    return res


def decode_column(number):
    res = []
    mask = 1
    for i in xrange(number.bit_length()):
        if number & mask:
            res.append(i)
        mask <<= 1
    return res

def ones(number):
    mask = 1
    res = 0
    for i in xrange(number.bit_length()):
        if number & mask:
            res += 1
        mask <<= 1
    return res


def bicluzt(mat, min_rows, min_cols):

    bitmat = bit_encode(mat)
    clusters = {}
    seen = set()
    for i in xrange(len(bitmat)):
        if bitmat[i] not in seen and ones(bitmat[i]) >= min_cols:
            i_clust = {bitmat[i] :  set([i])}

            for j in xrange(i+1, len(bitmat)):
                for seed in i_clust.keys():
                    j_seed = bitmat[j] & seed

                    if j_seed not in seen and ones(j_seed) >= min_cols:
                        i_clust.setdefault(j_seed, set())
                        i_clust[j_seed].update(i_clust[seed])
                        i_clust[j_seed].add(j)
            print len(i_clust)
            for seed in i_clust:
                seen.add(seed)
                if len(i_clust[seed]) >= min_rows:
                    clusters[seed] = i_clust[seed]

    return clusters

def BiBit(mat, min_rows, min_cols):

    bitmat = bit_encode(mat)
    clusters = {}
    for i in xrange(len(bitmat)):
        for j in xrange(i+1, len(bitmat)):
            p_ij = bitmat[i] & bitmat[j]
            if p_ij not in clusters and ones(p_ij) >= min_cols:
                cluster = set([i,j])
                for r in xrange(len(bitmat)):
                    if r == i or r == j:
                        continue
                    if bitmat[r] & p_ij == p_ij:
                        cluster.add(r)
                if len(cluster) >= min_rows:
                    clusters[p_ij] = cluster

    return clusters


def gwas(cases, controls, MAFs, alpha):
    res = []
    N = len(cases)
    M = len(cases[0])
    from scipy.stats import norm
    Scrit = norm().isf(alpha/(2*M))
    for snp in xrange(M):
        cases_freq = 0
        controls_freq = 0
        for ind in xrange(N):
            if cases[ind][snp]:
                cases_freq += 1
            if controls[ind][snp]:
                controls_freq += 1
        Sa = (float(cases_freq - controls_freq)/N)/(math.sqrt((2./N)*MAFs[snp]*(1-MAFs[snp])))
        if abs(Sa) > Scrit:
            res.append(snp)

    return res





test_mat = \
   [[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0],
    [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1],
    [0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1],
    [0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1],
    [0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0],
    [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
    [0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
    [0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0],
    [1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0],
    [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0]]


#[1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0],
#[1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1],
#[0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0],
#[0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
#[0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0],


#implanted: [([2, 4, 8, 13, 16], [3, 4, 7, 13, 16])]

SNPs = 500000
CASES = 1000
CONTROLS = 1000
BICLUSTERS = 10
BI_MAX_SNPs = 5
BI_MAX_INDs = 70

if __name__ == '__main__':
    case_m = matrix(CASES, SNPs)
    cont_m = matrix(CONTROLS, SNPs)
    implanted = []
    for bc in xrange(BICLUSTERS):
        bc_i = BI_MAX_INDs #random.randint(BI_MAX_INDs - 2, BI_MAX_INDs)
        bc_s = BI_MAX_SNPs #random.randint(BI_MAX_SNPs - 2, BI_MAX_SNPs)

        case_inds = sorted(random.sample(xrange(CASES), bc_i))
        case_snps = sorted(random.sample(xrange(SNPs), bc_s))
        implanted.append((case_snps, case_inds))

        for i in case_inds:
            for j in case_snps:
                case_m[i][j] = 1

    # add some noise
    MAFs = [0.1 for i in xrange(SNPs)]
    for j, maf in enumerate(MAFs):
        for i in xrange(CASES):
            if random.random() < maf:
                case_m[i][j] = 1
                cont_m[i][j] = 1

#    print pformat(case_m)
#    biclusters = BiBit(case_m, 5, 5)
#    biclusters = bicluzt(case_m, 5, 5)
#    print 'found:', pformat([(decode_column(key), individuals) for key, individuals in biclusters.iteritems()])
    print 'implanted:', pformat([imp[0] for imp in implanted])
    print 'gwas:', pformat(gwas(case_m, cont_m, MAFs, 0.05))
    elapsed('done')
#       control_inds = random.sample(xrange(CONTROLS))

