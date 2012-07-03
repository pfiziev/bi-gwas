import json
import os
import sys

__author__ = 'pf'


if __name__ == '__main__':
#    input_fname = '../random_GWAS.json'
#    input_fname = '../CEU_GWAS.json'

    input_fname = sys.argv[1]
    print input_fname

    data = json.load(open(input_fname))
    outf = open(os.path.split(input_fname)[1]+'.debi', 'w')

#    outf.write('ind_ids\t'+'\t'.join('snp_%d' % i for i in xrange(len(data['cases'][0])))+'\n')
#    outf.write('\n'.join('ind_%d\t' % i + '\t'.join(str(0 if v == 0 else 1) for v in row) for i, row in enumerate(data['cases'])))

    total_snps = len(data['cases'][0])
    total_inds = len(data['cases'])

    outf.write('snp_ids\t'+'\t'.join('ind_%d' % i for i in xrange(total_inds))+'\n')
    for snp_id in xrange(total_snps):
        outf.write('snp_%d\t' % snp_id + '\t'.join(str(0 if data['cases'][ind_id][snp_id] == 0 else 1) for ind_id in xrange(total_inds)) + '\n')


    outf.close()