import json
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


def binarize_matrix(mat):
    return [[int(bool(v)) for v in row] for row in mat]


if __name__ == '__main__':

#    data = json.load(open('CEU_GWAS.json'))
    data = json.load(open('random_GWAS.json'))

    cases = binarize_matrix(data['cases'])
    controls = binarize_matrix(data['controls'])
    implanted_biclusters = data['implanted_biclusters']

    MAFs = [(1 + sum(controls[ind][snp_id] for ind in xrange(len(controls))))/float(2*len(controls))
                for snp_id in xrange(len(controls[0]))]



#    print pformat(case_m)
#    biclusters = BiBit(cases, 5, 5)
    biclusters = bicluzt(cases, 5, 5)

    print 'found:', pformat([(decode_column(key), individuals) for key, individuals in biclusters.iteritems()])

    print 'implanted:', pformat([imp[1] for imp in implanted_biclusters])
    print 'gwas:', pformat(gwas(cases, controls, MAFs, 0.05))

    elapsed('done')
#       control_inds = random.sample(xrange(CONTROLS))

