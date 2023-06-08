import RNA

# sizes = [50, 100, 200, 500, 1000]
sizes = [50, 100, 200, 500, 1000]
nb = 5000
rhos = [1, 2, 5, 10, 50, 100, 500, 1021, 5000]


with open('ss_count_detail.csv', 'w') as g:
    print('Size', 'Salt', '# base pair', '# hairpin', '# multi loop', sep=',', file=g)
    for size in sizes:
        for rho in rhos:
            bps_counts = []
            hps_counts = []
            multis_counts = []
            with open('ss/ss_{}_{}.txt'.format(size, rho)) as f:
                for line in f.readlines():
                    ss = line.strip()
                    bps_counts.append(ss.count('('))
                    tree = RNA.db_to_tree_string(ss, RNA.STRUCTURE_TREE_SHAPIRO_SHORT)
                    hps_counts.append(tree.count('H'))
                    multis_counts.append(tree.count('M'))

            print(size, rho/1000, '{:.3f}'.format(sum(bps_counts)/nb), '{:.3f}'.format(sum(hps_counts)/nb), '{:.3f}'.format(sum(multis_counts)/nb), sep=',', file=g)
