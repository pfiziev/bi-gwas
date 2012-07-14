from itertools import izip
import json
import os
import random
import shlex
import shutil
import subprocess
import sys
import datetime
import math
#import numpy as np
import cPickle as pickle

__author__ = 'pf'


def now():
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

def bit_encode(mat):
    res = [0]*len(mat)
    for i, row in enumerate(mat):
        mask = 1
        for cell in row:
            if cell:
                res[i] += mask
            mask <<= 1

    return res

#WORD_SIZE = 24
#ones_counts = [-1] * 2**WORD_SIZE
#
#def ones(number):
##    res = ones.cache.get(number, None)
#
##    if not res:
#    mask = 1
#    res = 0
#    for i in xrange(number.bit_length()):
#        if number & mask:
#            res += 1
#        mask <<= 1
##        ones.cache[number] = res
#
#    return res
#ones.cache = {}
#


def ones(num):
    c = 0
    while num:
        num &= num - 1
        c+= 1
    return c


def find_potential_snps(cases, controls):
    CHUNK_SIZE = 20000
    potential_snps = set()

    print 'case snps:', len(cases)
    total_controls = float(len(controls[0]))
    print 'total controls:', total_controls
    potential_snps_genotypes = []

    for chunk_no in xrange(1+(len(cases) - 1)/CHUNK_SIZE):
        base_snp_id = chunk_no*CHUNK_SIZE
        last_snp_id = (chunk_no+1)*CHUNK_SIZE
        local_real = real_case_snps & set(xrange(base_snp_id, last_snp_id))

        print now(), 'processing chunk:', chunk_no, 'real snps:', len(local_real)
        print 'local real:', sorted(local_real)

        cases_bitmat = bit_encode(cases[base_snp_id: last_snp_id])

        cases_bitmat.extend(bit_encode(potential_snps_genotypes))

        print now(), 'bit encoding done'

        control_means = [sum(snp)/total_controls for snp in controls[base_snp_id: last_snp_id]]
        control_denominators = [math.sqrt(total_controls*p*(1-p)) for p in control_means]
        print now(), 'LD precalculations done'


        pairs = {}
        for i in xrange(len(cases_bitmat)):
            if i % 200 == 0: print now(), base_snp_id + i, len(pairs)

            if i < CHUNK_SIZE and (control_denominators[i] == 0 or ones(cases_bitmat[i]) < MIN_SUBJECTS):
                continue

            for j in xrange(i+1, len(cases_bitmat)):
                if (i < CHUNK_SIZE and j < CHUNK_SIZE) and \
                  (control_denominators[j] == 0 or
                   sum((x - control_means[i])*(y - control_means[j])
                            for x, y in izip(controls[i],controls[j]))/(control_denominators[i]*control_denominators[j]) >= 0.05):
                    continue

                if i >= CHUNK_SIZE and j >= CHUNK_SIZE:
                    continue

                pair = set()
                if i < CHUNK_SIZE:
                    pair.add(base_snp_id + i)
                if j < CHUNK_SIZE:
                    pair.add(base_snp_id + j)


                p_ij = cases_bitmat[i] & cases_bitmat[j]
#                x = ones(p_ij) >= MIN_SUBJECTS
#                continue
                if p_ij in pairs:
                    pairs[p_ij] |= pair

                elif ones(p_ij) >= MIN_SUBJECTS:
                    pairs[p_ij] = pair

        print now(), 'pairs found:', len(pairs)
        print 'local real:', sorted(local_real)
#        print pairs

        pair_keys = pairs.keys()

        for i in xrange(len(pair_keys)):
            for j in xrange(i+1, len(pair_keys)):
                if ones(pair_keys[i] & pair_keys[j]) >= MIN_SUBJECTS:
                    potential_snps |= pairs[pair_keys[i]] | pairs[pair_keys[j]]


        potential_snps_genotypes = [cases[snp_id] for snp_id in potential_snps]

        print now(), 'potential snps:', len(potential_snps)
        print now(), 'among them real snps:', len(real_case_snps & potential_snps)
#        print 'potential:', sorted(potential_snps)
#        print 'among them real:', sorted(real_case_snps & potential_snps)

        print '+'*100, '\n'

    return potential_snps



def run_in_parallel(commands):

    commands = [(cmd[0], open(cmd[1], 'w')) if type(cmd) is tuple else (cmd, None) for cmd in commands]

    print now(), 'Starting:\n' + '\n'.join([cmd for cmd, stdout in commands])
    for i, proc in enumerate([subprocess.Popen(args = shlex.split(cmd), stdout = stdout) for cmd, stdout in commands]):
        return_code = proc.wait()
        print now(), 'Finished: ' + commands[i][0]
        if return_code != 0:
            print now(), '%s \nexited with an error code: %d. Please, check the log files.' % (commands[i], return_code)
            exit(1)
    for _, stdout in commands:
        if stdout is not None:
            stdout.close()

if __name__ == '__main__':

    start_time = datetime.datetime.now()

    CHUNK_SIZE = 5000
    MIN_SUBJECTS = 90
#    input_fname = sys.argv[1]

#    input_fname = '../random_GWAS_300k.json'
#    input_fname = '../CEU_GWAS.json'

    input_fname = '../SIMLD/CEU_300k_chunked/CEU_300k.pickle'

    print now(), input_fname

#    data = json.load(open(input_fname))
    data = pickle.load(open(input_fname))

    tmp_dir = os.path.split(input_fname)[1] + '-tmp'
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.mkdir(tmp_dir)


    #    outf.write('ind_ids\t'+'\t'.join('snp_%d' % i for i in xrange(len(data['cases'][0])))+'\n')
    #    outf.write('\n'.join('ind_%d\t' % i + '\t'.join(str(0 if v == 0 else 1) for v in row) for i, row in enumerate(data['cases'])))

    total_snps = len(data['cases'])
    total_inds = len(data['cases'][0])
    print 'total snps:', total_snps
    print 'total individuals:', total_inds

    real_case_snps, real_case_inds = data['implanted_biclusters'][0]

    print 'real case inds:', len(real_case_inds)
    print 'real case snps:', len(real_case_snps)

    real_case_snps = set(real_case_snps)
    real_case_inds = set(real_case_inds)

    potential_snps = find_potential_snps(data['cases'], data['controls'])
    final_snps = set()

    print now(), 'initial potential snps:', len(potential_snps)
    print 'real snps among them:', len(potential_snps & real_case_snps)

    iteration = 0
#    total_iterations = 2*(1 + (total_snps - 1)/CHUNK_SIZE) + 1
    total_iterations = 1

    while iteration < total_iterations:
        print now(), 'iteration:', iteration, 'out of', total_iterations

#        snps_to_test = set(random.sample(xrange(total_snps), CHUNK_SIZE)) | potential_snps

        start_snp = (iteration % (total_snps/CHUNK_SIZE))*CHUNK_SIZE
        end_snp = min(total_snps, ((iteration%(total_snps/CHUNK_SIZE)) + 1)*CHUNK_SIZE)

        if iteration == total_iterations - 1:
#            chunk_snps = final_snps
            chunk_snps = potential_snps
            start_snp, end_snp = ['final', 'final']
        else:
            chunk_snps = set(xrange(start_snp, end_snp))


        print 'start:', start_snp, 'end:', end_snp
        print 'real snps:', len(chunk_snps & real_case_snps)
        print 'not in potential:', len((chunk_snps & real_case_snps) - potential_snps)
        print 'not in final:', len((chunk_snps & real_case_snps) - final_snps)

        snps_to_test = chunk_snps | potential_snps

        chunk_fname = os.path.join(tmp_dir, 'chunk.debi')

        outf = open(chunk_fname, 'w')
        outf.write('snp_ids\t'+'\t'.join('ind_%d' % i for i in xrange(total_inds))+'\n')

        for snp_id in snps_to_test:
            outf.write('snp_%d\t' % snp_id + '\t'.join(str(data['cases'][snp_id][ind_id]) for ind_id in xrange(total_inds)) + '\n')

        outf.close()

        run_in_parallel(['src/debi %s %s -o0.5 -pu -s%d' % (chunk_fname, tmp_dir, MIN_SUBJECTS) ])


        # check results
        tp = 0
        success = []

        print now(), 'checking results'

        lno = 0

        ind_cut = None


        for l in open(chunk_fname + '.biclusters'):
            lno += 1
            if lno == 2:
                snps = [int(v.split('_')[1]) for v in l.split()]

            elif lno == 3:
                inds = [int(v.split('_')[1]) for v in l.split()]

                ind_cut = set(inds) if ind_cut is None else ind_cut & set(inds)

                potential_snps = potential_snps | set(snps)
                final_snps = final_snps | set(snps)

                ind_overlap = len(set(inds) & real_case_inds) / float(len(real_case_inds))
                snp_overlap = len(set(snps) & real_case_snps) / float(len(real_case_snps))

                success.append((snp_overlap, ind_overlap, len(snps), len(inds)))

                lno = 0

        print 'Overlaps:', sorted(success, reverse=True)
        print 'potential snps:' , len(potential_snps), len(potential_snps & real_case_snps)
        print 'final snps:' , len(final_snps), len(final_snps & real_case_snps)
        if ind_cut:
            print 'inds:' , len(ind_cut), 'cut:', len(ind_cut & real_case_inds), len(real_case_inds)

        iteration += 1

    print 'elapsed:', datetime.datetime.now() - start_time



