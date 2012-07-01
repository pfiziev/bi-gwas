import json
import sys

if __name__ == '__main__':
    data = json.load(open(sys.argv[1]))
    tp = 0
    success = []
    lno = 0
    for l in open(sys.argv[2]):
        lno += 1
        if lno == 2:
            inds = [int(v.split('_')[1]) for v in l.split()]
        elif lno == 3:
            snps = [int(v.split('_')[1]) for v in l.split()]
#            print inds
#            print snps
#            print '+'*50
            for case_inds, case_snps in data['implanted_biclusters']:
#                print case_inds
#                print case_snps
                ind_overlap = len(set(inds) & set(case_inds)) / float(len(case_inds))
                snp_overlap = len(set(snps) & set(case_snps)) / float(len(case_snps))

                if snp_overlap > 95./100 and ind_overlap > 95./100:
                    tp += 1
                success.append((ind_overlap, snp_overlap))

            lno = 0
    print 'True positives:', tp
    print 'Overlaps:', success