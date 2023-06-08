from multiprocessing import Pool

import RNA

sizes = [50, 100, 200, 500, 1000]
nb = 5000
rhos = [1, 2, 5, 10, 50, 100, 500, 1021, 5000]


for size in sizes:
    print(size)
    seqLst = []
    with open('ss/seq_{}.txt'.format(size), 'w') as g:
        for _ in range(nb):
            w = RNA.random_string(size, 'ACGU')
            print(w, file=g)
            seqLst.append(w)
    for rho in rhos:
        print(rho)
        md = RNA.md(salt=rho/1000)

        def fold(w):
            fc = RNA.fold_compound(w, md)
            mfe, _ = fc.mfe()
            return mfe

        with open('ss/ss_{}_{}.txt'.format(size, rho), 'w') as f:
            with Pool(processes=8) as pool:
                for mfe in pool.map(fold, seqLst):
                    print(mfe, file=f)
