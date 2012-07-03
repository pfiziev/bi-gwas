import json
import sys

if __name__ == '__main__':
    data = json.load(open(sys.argv[1]))
    tp = 0
    success = []
    lno = 0

    #snp_cut = None
    #ind_union = set()
    
    snp_union = set()
    ind_cut = None
    
    for l in open(sys.argv[2]):
        lno += 1
        if lno == 2:
#            inds = [int(v.split('_')[1]) for v in l.split()]
            snps = [int(v.split('_')[1]) for v in l.split()]

        elif lno == 3:
#            snps = [int(v.split('_')[1]) for v in l.split()]
            inds = [int(v.split('_')[1]) for v in l.split()]

#            print inds
#            print snps
#            print '+'*50

#            snp_cut = set(snps) if snp_cut is None else snp_cut & set(snps)
#            ind_union = ind_union | set(inds)

            ind_cut = set(inds) if ind_cut is None else ind_cut & set(inds)
            snp_union = snp_union | set(snps)


            for case_inds, case_snps in data['implanted_biclusters']:
#                print case_inds
#                print case_snps
                ind_overlap = len(set(inds) & set(case_inds)) / float(len(case_inds))
                snp_overlap = len(set(snps) & set(case_snps)) / float(len(case_snps))

                if snp_overlap > 95./100:
                    tp += 1
                success.append((snp_overlap, ind_overlap, len(snps), len(inds)))

            lno = 0
    print 'True positives:', tp
    print 'Overlaps:', sorted(success, reverse=True)

#    print 'snps:' , len(snp_cut)
#    print 'inds:' , len(ind_union), 'cut:', len(ind_union & set(data['implanted_biclusters'][0][0])), len(data['implanted_biclusters'][0][0])

    print 'snps:' , len(snp_union), len(snp_union & set(data['implanted_biclusters'][0][1]))
    print 'inds:' , len(ind_cut), 'cut:', len(ind_cut & set(data['implanted_biclusters'][0][0])), len(data['implanted_biclusters'][0][0])

