import copy
import json
import sys
import datetime
import math

from pprint import pformat
import random
from itertools import *

import cPickle as pickle

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)

def now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


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

def ones(num):
    c = 0
    while num:
        num &= num - 1
        c+= 1
    return c

def bicluzt(mat, min_rows, min_cols):

    bitmat = bit_encode(mat)
    clusters = {}
    seen = set()
    for i in xrange(len(bitmat)):
        if bitmat[i] not in seen and ones(bitmat[i]) >= min_cols:
            elapsed('i: %d' % i)
            i_clust = {bitmat[i] :  set([i])}

            for j in xrange(i+1, len(bitmat)):
                for seed in i_clust.keys():
                    j_seed = bitmat[j] & seed

                    if j_seed not in seen and ones(j_seed) >= min_cols:
                        i_clust.setdefault(j_seed, set())
                        i_clust[j_seed].update(i_clust[seed])
                        i_clust[j_seed].add(j)

#                elapsed('j:%d\ti_clust length: %d' % (j, len(i_clust)))
            elapsed('i:%d\ti_clust length: %d' % (i, len(i_clust)))
            for seed in i_clust:
                seen.add(seed)
                if len(i_clust[seed]) >= min_rows:
                    clusters[seed] = i_clust[seed]

    return clusters

def BiBit(bitmat, min_rows, min_cols):

#    bitmat = bit_encode(mat)
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

        elapsed('i: %d' % i)

    return clusters


def find_potential_snps(cases_bitmat, controls_bitmat, potential_snps, potential_snp_genotypes, local_real, chunk_no, snp_mapping = None):
    total_controls = float(TOTAL_CONTROLS)
    total_snps = len(cases_bitmat)
    chunk_size = len(controls_bitmat)

    print now(), 'real snps:', len(local_real)


    print 'local real:', sorted(local_real)



#    cases_bitmat = bit_encode(cases)

    if len(potential_snps) > 0:
        pot_snp_ids, pot_snp_genotypes = zip(*potential_snp_genotypes.iteritems())
        cases_bitmat.extend(pot_snp_genotypes)
        total_snps += len(pot_snp_ids)

        if snp_mapping is not None:
            snp_mapping.extend(pot_snp_ids)



#    controls_bitmat = bit_encode(controls)


    print now(), 'bit encoding done'

    control_mafs = [ones(snp)/total_controls for snp in controls_bitmat]
    control_denominators = [math.sqrt(p*(1-p)) for p in control_mafs]

    print now(), 'LD precalculations done'


    pairs = {}
    for i in xrange(chunk_size):
        if i % 200 == 0:
            print now(), i, snp_mapping[i], len(pairs)
            if len(pairs) > chunk_size:
                print "Giving up!"
                break

        if control_denominators[i] == 0 or ones(cases_bitmat[i]) < MIN_SUBJECTS:
            continue

        for j in xrange(i+1, total_snps):
            if j < chunk_size and \
               (control_denominators[j] == 0 or

                (ones(controls_bitmat[i] & controls_bitmat[j])/total_controls - control_mafs[i]*control_mafs[j])/
                (control_denominators[i]*control_denominators[j]) >= MAX_LD):

                continue


            p_ij = cases_bitmat[i] & cases_bitmat[j]
            #                x = ones(p_ij) >= MIN_SUBJECTS
            #                continue
            if p_ij in pairs or ones(p_ij) >= MIN_SUBJECTS:

                pair = { snp_mapping[i],
                         snp_mapping[j] }


                if p_ij in pairs:

                    pairs[p_ij] |= pair
                else:
                    pairs[p_ij] = pair

    print now(), 'pairs found:', len(pairs)
    print 'local real:', sorted(local_real)
    print 'pairs:', pairs.values() if len(pairs) < 100 else '> 100'

    pair_keys = pairs.keys()

    for i in xrange(len(pair_keys)):
        if len(pairs[pair_keys[i]]) >= 3:
            potential_snps |= pairs[pair_keys[i]]

        for j in xrange(i+1, len(pair_keys)):
            if ones(pair_keys[i] & pair_keys[j]) >= MIN_SUBJECTS:
                potential_snps |= pairs[pair_keys[i]] | pairs[pair_keys[j]]


    print now(), 'potential snps:', len(potential_snps)
    print now(), 'local not found:', len(local_real - potential_snps)
    #        print 'potential:', sorted(potential_snps)
    #        print 'among them real:', sorted(real_case_snps & potential_snps)

    print 'potential snps list:', potential_snps if len(potential_snps) < 100 else '>100'
    return potential_snps



def preprocess(cases, controls, real_case_snps):
    CHUNK_SIZE = 10000
    potential_snps = set()
    potential_snp_genotypes = {}

    total_snps = len(cases)

    print 'chunk size:', CHUNK_SIZE
    print 'total snps:', total_snps
    print 'causal snps:', len(real_case_snps)
    total_controls = float(TOTAL_CONTROLS)
    print 'total controls:', total_controls
    print 'total cases:', TOTAL_CASES

    chunks_to_sample = how_many_chunks(total_snps , CHUNK_SIZE, MIN_SNPS, 0.95)
    print now(), 'chunks to sample:', chunks_to_sample

    for chunk_no in xrange(chunks_to_sample):
        print now(), 'sampling random chunk:', chunk_no, 'out of', chunks_to_sample

        snp_mapping = sorted(random.sample(xrange(total_snps), CHUNK_SIZE))

        sample_cases = [cases[snp_id] for snp_id in snp_mapping]
        sample_controls= [controls[snp_id] for snp_id in snp_mapping]

        potential_snp_genotypes = dict((snp_id, cases[snp_id]) for snp_id in potential_snps)

        find_potential_snps(sample_cases,
                            sample_controls,
                            potential_snps,
                            potential_snp_genotypes,
                            real_case_snps & set(snp_mapping),
                            chunk_no,
                            snp_mapping = snp_mapping )

        print now(), 'potential snps:',len(potential_snps)#, sorted(potential_snps)
        print now(), 'among them real snps:', len(real_case_snps & potential_snps)
        print '+'*100, '\n'


    chunks_to_process = 1+(total_snps - 1)/CHUNK_SIZE
    for chunk_no in xrange(chunks_to_process):

        base_snp_id = chunk_no*CHUNK_SIZE
        last_snp_id = (chunk_no+1)*CHUNK_SIZE
        snp_mapping = range(base_snp_id, last_snp_id)
        local_real  = real_case_snps & set(snp_mapping)

        potential_snp_genotypes = dict((snp_id, cases[snp_id]) for snp_id in potential_snps)
        print now(), 'processing chunk:', chunk_no, 'out of', chunks_to_process

        find_potential_snps(  cases[base_snp_id:last_snp_id],
                              controls[base_snp_id:last_snp_id],
                              potential_snps,
                              potential_snp_genotypes,
                              local_real,
                              chunk_no,
                              snp_mapping = snp_mapping)


        print now(), 'potential snps:', len(potential_snps), sorted(potential_snps)
        print now(), 'among them real snps:', len(real_case_snps & potential_snps)
        print '+'*100, '\n'


    print now(), 'not found:', len(real_case_snps - potential_snps), real_case_snps - potential_snps


    return potential_snps

#
#
#def _preprocess(cases, controls, real_case_snps):
#    CHUNK_SIZE = 10000
#    potential_snps = set()
#    total_snps = len(cases)
#
#    print 'chunk size:', CHUNK_SIZE
#    print 'total snps:', total_snps
#    print 'causal snps:', len(real_case_snps)
#    total_controls = float(len(controls[0]))
#    print 'total controls:', total_controls
#    potential_snps_genotypes = []
#
#    chunks_to_sample = how_many_chunks(total_snps , CHUNK_SIZE, len(real_case_snps), 0.95)
#    print now(), 'chunks to sample:', chunks_to_sample
#    for chunk_no in xrange(chunks_to_sample):
#        print now(), 'sampling random chunk:', chunk_no
#
#        sample_snps = random.sample(xrange(total_snps), CHUNK_SIZE)
#
#        sample_cases = [cases[snp_id] for snp_id in sample_snps]
#        sample_cases.extend([cases[snp_id] for snp_id in potential_snps])
#
#        sample_controls= [controls[snp_id] for snp_id in sample_snps]
#
#        potential_snps |= find_potential_snps(sample_cases, sample_controls)
#
#
#    for chunk_no in xrange(1+(len(cases) - 1)/CHUNK_SIZE):
#        base_snp_id = chunk_no*CHUNK_SIZE
#        last_snp_id = (chunk_no+1)*CHUNK_SIZE
#        local_real = real_case_snps & set(xrange(base_snp_id, last_snp_id))
#
#        print now(), 'processing chunk:', chunk_no, 'real snps:', len(local_real)
#        print 'local real:', sorted(local_real)
#        print 'new local real:' , len(local_real - potential_snps), local_real - potential_snps
#
#        cases_bitmat = bit_encode(cases[base_snp_id: last_snp_id])
#
#        cases_bitmat.extend(bit_encode(potential_snps_genotypes))
#
#        print now(), 'bit encoding done'
#
##        control_means = [sum(snp)/total_controls for snp in controls[base_snp_id: last_snp_id]]
##        control_denominators = [math.sqrt(total_controls*p*(1-p)) for p in control_means]
#
#        controls_bitmat = bit_encode(controls[base_snp_id: last_snp_id])
#
#        control_mafs = [ones(snp)/total_controls for snp in controls_bitmat]
#        control_denominators = [math.sqrt(p*(1-p)) for p in control_mafs]
#
#        print now(), 'LD precalculations done'
#
#
#        pairs = {}
#        for i in xrange(len(cases_bitmat)):
#            if i % 200 == 0: print now(), base_snp_id + i, len(pairs)
#
#            if i < CHUNK_SIZE and (control_denominators[i] == 0 or ones(cases_bitmat[i]) < MIN_SUBJECTS):
#                continue
#
#            for j in xrange(i+1, len(cases_bitmat)):
#                if (i < CHUNK_SIZE and j < CHUNK_SIZE) and\
#                   (control_denominators[j] == 0 or
#
#                   (ones(controls_bitmat[i] & controls_bitmat[j])/total_controls - control_mafs[i]*control_mafs[j])/
#                   (control_denominators[i]*control_denominators[j]) >= 0.2
#
##                    sum((x - control_means[i])*(y - control_means[j])
##                        for x, y in izip(controls[i],controls[j]))/(control_denominators[i]*control_denominators[j]) >= 0.05
#
#                       ):
#                    continue
#
#                if i >= CHUNK_SIZE and j >= CHUNK_SIZE:
#                    continue
#
#                pair = set()
#                if i < CHUNK_SIZE:
#                    pair.add(base_snp_id + i)
#                if j < CHUNK_SIZE:
#                    pair.add(base_snp_id + j)
#
#
#                p_ij = cases_bitmat[i] & cases_bitmat[j]
#                #                x = ones(p_ij) >= MIN_SUBJECTS
#                #                continue
#                if p_ij in pairs:
#                    pairs[p_ij] |= pair
#
#                elif ones(p_ij) >= MIN_SUBJECTS:
#                    pairs[p_ij] = pair
#
#        print now(), 'pairs found:', len(pairs)
#        print 'local real:', sorted(local_real)
#        #        print pairs
#
#        pair_keys = pairs.keys()
#
#        for i in xrange(len(pair_keys)):
#            for j in xrange(i+1, len(pair_keys)):
#                if ones(pair_keys[i] & pair_keys[j]) >= MIN_SUBJECTS:
#                    potential_snps |= pairs[pair_keys[i]] | pairs[pair_keys[j]]
#
#
#        potential_snps_genotypes = [cases[snp_id] for snp_id in potential_snps]
#
#        print now(), 'potential snps:', len(potential_snps)
#        print now(), 'among them real snps:', len(real_case_snps & potential_snps)
#        print now(), 'local not found:', len(local_real - potential_snps)
#        #        print 'potential:', sorted(potential_snps)
#        #        print 'among them real:', sorted(real_case_snps & potential_snps)
#
#        print '+'*100, '\n'
#
#
#
#    return potential_snps



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



def pval(snps, cases, controls):
    return 0
#    bi_controls = ones(reduce(lambda x,y: x & y, bit_encode([controls[snp_id] for snp_id in snps])))
#    bi_cases    = ones(reduce(lambda x,y: x & y, bit_encode([cases[snp_id] for snp_id in snps])))
#    return 'controls: %d' % bi_controls



def how_many_chunks(total_snps, chunk_size, min_causal_snps, alpha):
    """
    Calculates how many chunks we need to consider before being at least alpha % confident
    that at least one chunk contained minimum 3 causal SNPs
    """
    binom_cdf_2 = lambda n, p:  (1-p)**n + n*((1-p)**(n-1))*p + (n*(n-1)/2)*((1-p)**(n-2))*(p**2)
    total_chunks = total_snps/chunk_size
    P_c = (1 - binom_cdf_2(min_causal_snps, 1./total_chunks))
    return int(math.ceil(math.log(1-alpha, 1 - P_c)))



if __name__ == '__main__':


    #    data = json.load(open('CEU_GWAS.json'))
#    data = json.load(open('random_GWAS2.json'))
#    cases = binarize_matrix(data['cases'])
#    controls = binarize_matrix(data['controls'])
#    implanted_biclusters = data['implanted_biclusters']
#    input_file = 'SIMLD/CEU_300k_10k_chunked/CEU_100k_30_SNPs_by_100_INDS.pickle'


#    input_file = 'SIMLD/CEU_300k_10k_chunked/CEU_100k_10_SNPs_by_100_INDS.pickle'

#    input_file = 'SIMLD/CEU_300k_10k_chunked/CEU_100k_30_SNPs_by_80_INDS.pickle'

#    input_file = 'SIMLD/CEU_300k_10k_chunked/CEU_100k_30_SNPs_by_70_INDS.pickle'

#    input_file = 'SIMLD/CEU_300k_10k_chunked/CEU_300k_100_SNPs_by_100_INDS.pickle'
    input_file = 'SIMLD/CEU_300k_10k_chunked/CEU_300k_100_SNPs_by_70_INDS.pickle'


#    input_file = 'SIMLD/CEU_300k_10k_chunked/CEU_100k_15_SNPs_by_100_INDS.pickle'
#    input_file = 'SIMLD/CEU_300k_10k_chunked/CEU_10k_30_SNPs_by_100_INDS.pickle'

    print 'input:', input_file

    cases, TOTAL_CASES, controls, TOTAL_CONTROLS, implanted_biclusters = pickle.load(open(input_file))

#    data = pickle.load(open(input_file))
#    cases = data['cases']
#    controls = data['controls']
#    implanted_biclusters = data['implanted_biclusters']

    MAX_LD = 0.2

    MIN_SUBJECTS = min(len(subj) for _, subj in implanted_biclusters)
    MIN_SNPS = min(len(bi_snps) for bi_snps, _ in implanted_biclusters)

    print 'MIN_SUBJECTS:', MIN_SUBJECTS
    print 'MIN_SNPs:', MIN_SNPS
    print 'MAX LD:', MAX_LD


    potential_snps = sorted(preprocess(cases, controls, set(implanted_biclusters[0][0])))


    elapsed('preprocessing')

    to_cluster = [cases[snp_id] for snp_id in potential_snps]

    to_cluster_implanted = copy.deepcopy(implanted_biclusters)

    for bi_snps, bi_inds in to_cluster_implanted:
        print 'original snp ids:', bi_snps

        for i in xrange(len(bi_snps)):
            bi_snps[i] = potential_snps.index(bi_snps[i]) if bi_snps[i] in potential_snps else -1

    pickle.dump({'to_cluster' : to_cluster, 'implanted' : to_cluster_implanted},
                open(input_file + '.preprocessed.pickle', 'w'))

    print 'final implanted_biclusters:', to_cluster_implanted
#    MAFs = [(1 + sum(controls[ind][snp_id] for ind in xrange(len(controls))))/float(2*len(controls))
#                for snp_id in xrange(len(controls[0]))]


    elapsed('cases: %d\tcontrols: %d' % (len(cases), len(controls)))

#    biclusters = bicluzt(to_cluster, MIN_SNPS, MIN_SUBJECTS)

#    print pformat(case_m)
    biclusters = BiBit(to_cluster, MIN_SNPS, MIN_SUBJECTS)

    found_snps = [set(biclusters[key]) for key in biclusters]

#    print 'found:', pformat([(decode_column(key), individuals) for key, individuals in biclusters.iteritems()])
    print 'found:'
    for snps in found_snps:

        print 'pvalue:',pval(snps, cases, controls), 'count:', len(snps), 'overlap:', [len(snps & set(implanted)) for implanted, _ in to_cluster_implanted]

#    print 'implanted:', pformat([imp[0] for imp in to_cluster_implanted])
    elapsed('done biclustering!')

#    total_controls = len(controls[0])
#
#    MAFs = [sum(snp)/total_controls]
#    print 'gwas:', pformat(gwas(cases, controls, MAFs, 0.05))

    elapsed('done')
#       control_inds = random.sample(xrange(CONTROLS))



## code to calculate min number of trials
#
#total_snps = 500000
#chunk_size = 20000
#total_chunks= total_snps/chunk_size
#causal_snps = 20
#print 'total chunks:', total_chunks
## P_c = P(X >=3|n,p) = 1 - P(X <= 2| n,p) !!!
#P_c = (1-binom(causal_snps, 1./total_chunks).cdf(2))
#print 'P_c:', P_c
#print 'E(#chunks):', P_c * total_chunks
#tries = 70
#print 'alpha(%d):'%tries, 1-(1-P_c)**tries


## code to plot the allele distributions for cases and controls
#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
#
#import cPickle as pickle
#
#data = pickle.load(open('SIMLD/CEU_300k_10k_chunked/CEU_300k_100_SNPs_by_100_INDS.pickle'))
#
#maf_cases = [float(sum(snp))/len(snp) for snp in data['cases']]
#maf_controls = [float(sum(snp))/len(snp) for snp in data['controls']]
#
#
#fig = plt.figure()
#ax = fig.add_subplot(111)
#
## the histogram of the data
#n, bins, patches = ax.hist(maf_cases, 50, normed=1, facecolor='green', alpha=0.75)
#ax.set_title('Homozygous allele distribution in cases')
#
#fig = plt.figure()
#ax = fig.add_subplot(111)
#
## the histogram of the data
#n, bins, patches = ax.hist(maf_controls, 50, normed=1, facecolor='green', alpha=0.75)
#ax.set_title('Homozygous allele distribution in controls')


##code to calculate the LD structure:
#
#
#def ld(bitmat):
#    total_inds = 1000.
#    total_snps = len(bitmat)
#    mafs = [ones(snp)/total_inds for snp in bitmat]
#    #    print [m for m in  mafs if m > 1 or m < 0]
#    dens = [math.sqrt(p*(1-p)) for p in mafs]
#    ld_mat = zeros((total_snps, total_snps))
#    for i in xrange(total_snps):
#        ld_mat[i][i] = 1
#        if i > 0:
#            ld_mat[i][i-1] = -1
#        for j in xrange(i + 1, total_snps):
#            if dens[i] and dens[j]:
#                ld_mat[i][j] = (ones(bitmat[i] & bitmat[j])/total_inds - mafs[i]*mafs[j])/(dens[i]*dens[j])
#
#    return ld_mat
#
#
#from matplotlib.pylab import *
#
#
#def bit_encode(mat):
#    res = [0]*len(mat)
#    for i, row in enumerate(mat):
#        mask = 1
#        for cell in row:
#            if cell:
#                res[i] += mask
#            mask <<= 1
#
#    return res
#
#def ones(num):
#    c = 0
#    while num:
#        num &= num - 1
#        c+= 1
#    return c
#
#data = pickle.load(open('/home/pf/local_temp/bi-gwas/SIMLD/CEU_300k_10k_chunked/CEU_100k_30_SNPs_by_70_INDS.pickle'))
#cases = data['cases'][:200]
#controls = data['controls'][1000:1200]
#
#matshow(ld(bit_encode(cases)))
