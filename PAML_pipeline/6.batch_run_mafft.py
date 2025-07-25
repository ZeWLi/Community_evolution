'''

'''

import os
import sys

def batch_run_with_msa(dir,t):
    path = []
    for root,dirs,files in os.walk(dir):
        for file in files:
            path.append(file)

    os.chdir(dir)
    for i in path:
        os.system(f'linsi --thread {t} --maxiterate 1000 {i} > {i.strip(".fa")}.mafft.fa')

    os.chdir('..')
    return 0

def main():
    batch_run_with_msa(sys.argv[1],eval(sys.argv[2]))
    return 0

if __name__ =='__main__':
    if len(sys.argv) != 3:
        print('wrong argv')
        exit()

    main()