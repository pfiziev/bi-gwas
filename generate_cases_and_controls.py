"""
This script converts CEU.BEAM.txt (which is in BEAM format: first line is 0 or 1 if subject is case or control,
then one line per SNP and 0 if homozygous for the major allele, 1 if heterozygous, 2 if homozygous for the minor allele).

The output is two matrices: cases and controls
and a description of the implanted biclusters.

cases and controls a are transposed matrices:
Each row corresponds to a person. Each column corresponds to a SNP. Cells are 0, 1 or 2 (same meaning as in BEAM format)

SNPs in the biclusters are set to 2.
"""


import json
import os
import cPickle as pickle
import random
#from scipy.stats import pearsonr
from utils import *



def generate_random_cases_and_controls():

    SNPs = 300000
    CASES = 1000
    CONTROLS = 100
    BICLUSTERS = 1
    BI_MAX_SNPs = 101
    BI_MAX_INDs = 100

    case_m = matrix(SNPs, CASES)
    cont_m = matrix(SNPs, CONTROLS)

    implanted_biclusters = []

    # add some noise
    MAFs = [0.1] * SNPs

#    MAFs = [0.05 + (float(i)/(3*SNPs)) for i in xrange(SNPs)]

    for i, maf in enumerate(MAFs):
        for j in xrange(CASES):
            case_m[i][j] = (1 if random.random() < maf else 0) + (1 if random.random() < maf else 0)

        for j in xrange(CONTROLS):
            cont_m[i][j] = (1 if random.random() < maf else 0) + (1 if random.random() < maf else 0)


    for bc in xrange(BICLUSTERS):

        bc_i = BI_MAX_INDs #random.randint(BI_MAX_INDs - 2, BI_MAX_INDs)
        bc_s = BI_MAX_SNPs #random.randint(BI_MAX_SNPs - 2, BI_MAX_SNPs)

        case_inds = sorted(random.sample(xrange(CASES), bc_i))
        case_snps = sorted(random.sample(xrange(SNPs), bc_s))

        implanted_biclusters.append((case_snps, case_inds))
        for i in case_snps:
            for j in case_inds:
                case_m[i][j] = 2


    json.dump({'cases': case_m, 'controls': cont_m, 'implanted_biclusters' : implanted_biclusters}, open('random_GWAS_300k.json', 'w'))


def generate_from_BEAM():

    PPL_TO_TAKE = 1000

    CASE_CONTROL_RATIO = 0.5 # the ratio of cases vs controls

    BI_CASES = 100 # fraction of cases that are in one bicluster
    BI_SNPS = 100   # number of SNPs per bicluster
    BICLUSTERS = 1 # number of biclusters

    snp_file = open('SIMLD/CEU.BEAM.txt')

    # read out disease status. this information is irrelevant.
    _ = snp_file.readline()

    snps = [map(int, l.split()) for l in snp_file]
    snp_file.close()

    total_snps = len(snps)
    pop_size = len(snps[0])

    total_cases = int(CASE_CONTROL_RATIO * pop_size)

    # create cases and controls matrices and transpose them:
#    cases    = [list(row) for row in zip(*[snp[:total_cases] for snp in snps])][:PPL_TO_TAKE]
#    controls = [list(row) for row in zip(*[snp[total_cases:] for snp in snps])][:PPL_TO_TAKE]

    # don't transpose anything
    cases    = [snp[:total_cases] for snp in snps]
    controls = [snp[total_cases:] for snp in snps]

    total_cases = min(PPL_TO_TAKE, total_cases)

    implanted_biclusters = [[ random.sample(xrange(total_snps), BI_SNPS),
                              random.sample(xrange(total_cases), int(BI_CASES))]
                                        for _ in xrange(BICLUSTERS)]



    for bi_snps, bi_ppl in implanted_biclusters:
        print 'implanting - people:', len(bi_ppl), ', snps:', len(bi_snps)

        for snp_id in bi_snps:
            for person_id in bi_ppl:
                cases[snp_id][person_id] = 2

#    snp_freq = [sum(snp)/float(2*pop_size) for snp in snps]
#    cc = [[0 for i in xrange(total_snps)] for j in xrange(total_snps)]
#
#    for i in xrange(total_snps):
#        for j in xrange(total_snps):
#            cc[i][j] = pearsonr(snps[i], snps[j])[0]
#    print min(snp_freq), max(snp_freq)
#    print total_snps, pop_size

    json.dump({'cases': cases, 'controls': controls, 'implanted_biclusters' : implanted_biclusters}, open('CEU_GWAS.json', 'w'))



def generate_from_BEAM_chunks():
    HOMOZYGOUS = 1
    HETEROZYGOUS = 0

    CASE_CONTROL_RATIO = 0.5 # the ratio of cases vs controls

    BI_CASES = 100 # fraction of cases that are in one bicluster
    BI_SNPS = 30   # number of SNPs per bicluster
    BICLUSTERS = 1 # number of biclusters
    FILES_TO_TAKE = 20

    TOTAL_CASES = 1000
    TOTAL_INDIVIDUALS = 2000
    case_ids = random.sample(xrange(TOTAL_INDIVIDUALS), TOTAL_CASES)
    control_ids = [pid for pid in xrange(TOTAL_INDIVIDUALS) if pid not in case_ids]


    snp_dir = 'SIMLD/CEU_300k_chunked'

    cases = []
    controls = []

    for snp_fname in map(lambda f: os.path.join(snp_dir, f), sorted([f for f in os.listdir(snp_dir) if f.endswith('.txt')]))[:FILES_TO_TAKE]:
        print 'processing', snp_fname

        snp_file = open(snp_fname)

        # read out disease status to deterimin population size.
        pop_size = len(snp_file.readline().split())
        total_cases = int(CASE_CONTROL_RATIO * pop_size)

        for l in snp_file:
            snps = [0 if v == '0' else
                   (HETEROZYGOUS if v == '1' else
                   (HOMOZYGOUS if v == '2' else None)) for v in l.split()]

            cases.append([snps[pid] for pid in case_ids])
            controls.append([snps[pid] for pid in control_ids])


        snp_file.close()

    total_snps = len(cases)

    implanted_biclusters = [[ sorted(random.sample(xrange(total_snps), BI_SNPS)),
                              sorted(random.sample(xrange(total_cases), BI_CASES))]
                                    for _ in xrange(BICLUSTERS)]


    for bi_snps, bi_ppl in implanted_biclusters:
        print 'implanting - people:', len(bi_ppl), ', snps:', len(bi_snps)

        for snp_id in bi_snps:
            for person_id in bi_ppl:
                cases[snp_id][person_id] = HOMOZYGOUS



    pickle.dump({ 'cases': cases,
                  'controls': controls,
                  'implanted_biclusters' : implanted_biclusters},
                open('SIMLD/CEU_300k_chunked/CEU_300k.pickle', 'w'),
                pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
#    generate_from_BEAM()
    generate_from_BEAM_chunks()
#    generate_random_cases_and_controls()
