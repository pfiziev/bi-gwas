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
import random
from scipy.stats import pearsonr
from utils import *



def generate_random_cases_and_controls():

    SNPs = 30000
    CASES = 1000
    CONTROLS = 100
    BICLUSTERS = 1
    BI_MAX_SNPs = 100
    BI_MAX_INDs = 100

    case_m = matrix(CASES, SNPs)
    cont_m = matrix(CONTROLS, SNPs)
    implanted_biclusters = []

    # add some noise
    MAFs = [0.1] * SNPs

#    MAFs = [0.05 + (float(i)/(3*SNPs)) for i in xrange(SNPs)]

    for j, maf in enumerate(MAFs):
        for i in xrange(CASES):
            case_m[i][j] = (1 if random.random() < maf else 0) + (1 if random.random() < maf else 0)

        for i in xrange(CONTROLS):
            cont_m[i][j] = (1 if random.random() < maf else 0) + (1 if random.random() < maf else 0)


    for bc in xrange(BICLUSTERS):

        bc_i = BI_MAX_INDs #random.randint(BI_MAX_INDs - 2, BI_MAX_INDs)
        bc_s = BI_MAX_SNPs #random.randint(BI_MAX_SNPs - 2, BI_MAX_SNPs)

        case_inds = sorted(random.sample(xrange(CASES), bc_i))
        case_snps = sorted(random.sample(xrange(SNPs), bc_s))

        implanted_biclusters.append((case_inds, case_snps))

        for i in case_inds:
            for j in case_snps:
                case_m[i][j] = 2


    json.dump({'cases': case_m, 'controls': cont_m, 'implanted_biclusters' : implanted_biclusters}, open('random_GWAS.json', 'w'))


def generate_from_BEAM():
    PPL_TO_TAKE = 1000

    CASE_CONTROL_RATIO = 0.5 # the ratio of cases vs controls

    BI_CASES = 100 # fraction of cases that are in one bicluster
    BI_SNPS = 50   # number of SNPs per bicluster
    BICLUSTERS = 1 # number of biclusters

    snp_file = open('SIMLD/CEU.BEAM.txt')
    _ = snp_file.readline()
    snps = [map(int, l.split()) for l in snp_file]
    snp_file.close()
    total_snps = len(snps)
    pop_size = len(snps[0])

    total_cases = int(CASE_CONTROL_RATIO * pop_size)

    # create cases and controls matrices and transpose them:
    cases    = [list(row) for row in zip(*[snp[:total_cases] for snp in snps])][:PPL_TO_TAKE]
    controls = [list(row) for row in zip(*[snp[total_cases:] for snp in snps])][:PPL_TO_TAKE]

    total_cases = min(PPL_TO_TAKE, total_cases)

    implanted_biclusters = [[ random.sample(xrange(total_cases), int(BI_CASES)),
                    random.sample(xrange(total_snps), BI_SNPS)]
                        for _ in xrange(BICLUSTERS)]



    for bi_ppl, bi_snps in implanted_biclusters:
        print 'implanting - people:', len(bi_ppl), ', snps:', len(bi_snps)
        for person_id in bi_ppl:
            for snp_id in bi_snps:
                cases[person_id][snp_id] = 2

#    snp_freq = [sum(snp)/float(2*pop_size) for snp in snps]
#    cc = [[0 for i in xrange(total_snps)] for j in xrange(total_snps)]
#
#    for i in xrange(total_snps):
#        for j in xrange(total_snps):
#            cc[i][j] = pearsonr(snps[i], snps[j])[0]
#    print min(snp_freq), max(snp_freq)
#    print total_snps, pop_size

    json.dump({'cases': cases, 'controls': controls, 'implanted_biclusters' : implanted_biclusters}, open('CEU_GWAS.json', 'w'))



if __name__ == '__main__':
#    generate_from_BEAM()
    generate_random_cases_and_controls()
