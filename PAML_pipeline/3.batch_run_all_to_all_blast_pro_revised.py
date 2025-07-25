'''
Date:20230725
nohup python /data/lizw/script/paml_pipeline/batch_run_all_to_all_blast_pro.py -d 1.gene_test &
'''

import os
import argparse

def all_to_all_blast(t):
    os.mkdir('temp_database')
    os.chdir('temp_database')
    os.system('makeblastdb -in ../all_pro.trim33.fa -dbtype prot -out temp')
    os.chdir('..')
    os.system(f'blastp -query all_pro.trim33.fa -out blast.out -db temp_database/temp -evalue 1e-5 -num_threads {t} -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"')
    os.system('rm -r temp_database')
    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', help='#threads of blast',type=int,default=50)
    args = parser.parse_args()

    all_to_all_blast(args.t)