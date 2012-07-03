import json
import sys

if __name__ == '__main__':
    data = json.load(open(sys.argv[1]))
    tp = 0
    success = []

    for l in open(sys.argv[2]):
        bicl = json.loads(l)
        inds = bicl['inds']
        snps = bicl['snps']

        #            print inds
        #            print snps
        #            print '+'*50
        for case_inds, case_snps in data['implanted_biclusters']:
        #                print case_inds
        #                print case_snps
            ind_overlap = len(set(inds) & set(case_inds)) / float(len(case_inds))
            snp_overlap = len(set(snps) & set(case_snps)) / float(len(case_snps))

            if snp_overlap > 95./100:
                tp += 1
            success.append((snp_overlap, ind_overlap, len(snps), len(inds)))
    print 'True positives:', tp
    print 'Overlaps:', sorted(success, reverse=True)